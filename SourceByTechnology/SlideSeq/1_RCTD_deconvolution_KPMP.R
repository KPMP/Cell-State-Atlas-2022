library(ggplot2)
library(Matrix)
library(Giotto)
library(reshape2)
library(patchwork)
library(RCTD)
library(Seurat)
library(foreach)
source("./utils.R")

fitPixels_full <- function(RCTD){
  # Modified version of RCTD's "fitPixels" function that runs faster
  doublet_mode <- 'full'
  RCTD@config$doublet_mode <- doublet_mode
  set_likelihood_vars(RCTD@internal_vars$Q_mat, RCTD@internal_vars$X_vals)
  cell_type_info <- RCTD@cell_type_info$renorm
  
  if (doublet_mode == "full") {
    beads = t(as.matrix(RCTD@spatialRNA@counts[RCTD@internal_vars$gene_list_reg, 
                                               ]))
    results = RCTD:::decompose_batch(RCTD@spatialRNA@nUMI, cell_type_info[[1]], 
                                     beads, RCTD@internal_vars$gene_list_reg, constrain = F, 
                                     max_cores = RCTD@config$max_cores, MIN.CHANGE = RCTD@config$MIN_CHANGE_REG)
    weights <- sapply(results, function(x) x$weights)
    weights <- t(weights)
    rownames(weights) = colnames(RCTD@spatialRNA@counts)
    colnames(weights) = RCTD@cell_type_info$renorm[[2]]
    weights <- Matrix(weights, sparse = T)
    RCTD@results <- list(weights = weights)
  }
  return(RCTD)
}

# HuBMAP data comes in two batches that are in slightly different formats
datadirs <- c("/media/NAS2/Users/Kian_NAS2/SlideSeq_data/Broad_JamieMarshall/")
colorTable <- data.table::fread("./Cluster_Color_Table.tsv", select = 1:20)
python_path <- "/media/Home_Raid1_Voyager/kian/anaconda3/envs/giotto_200703/bin/python"
seurat.path <- "../../../kidney_references/Kidney_KPMP-Biopsy_Ref_10X-R_092020_Seurat_F.rda"

ex.lines <- list() # a list containing the tissue boundaries for some pucks

out.dir <- "./deconv_results"
if (!dir.exists(out.dir))
  dir.create(out.dir)

rctd.dir <- out.dir
if (!dir.exists(rctd.dir))
  dir.create(rctd.dir)

plot.dir <- paste0("./Plots")
if (!dir.exists(plot.dir))
  dir.create(plot.dir)

gioObjFolder <- paste0(rctd.dir, "/GiottoObjects")
if (!dir.exists(gioObjFolder))
  dir.create(gioObjFolder)


prefix.f <- "KPMP"
l1file <- paste0(rctd.dir, "/", prefix.f, "_l1.rds")
l2file <- paste0(rctd.dir, "/", prefix.f, "_l2.rds")

# # Log file
sink(paste0(out.dir, "/", Sys.time(), '_', prefix.f, ".log"))

# prepping puck names and addressed
pucks <- list.dirs(datadirs, full.names = F)
pucks <- pucks[nchar(pucks) > 0]
puckdirs <- paste(datadirs, pucks, sep = '/')
puckdirs <- as.list(puckdirs)
names(puckdirs) <- pucks
print(puckdirs)


# Reading all the pucks at once, making one Giotto object to contain them all, discarding the spatial locations
# Reading counts ####
umi.thresh <- 100
counts.all <- list()
for(puck in names(puckdirs)){
  print(paste('Reading', puck))
  counts.all[[puck]] <- readSlideSeqGiotto3(puck, datadirs, python_path=python_path)@raw_exprs
  counts.all[[puck]] <- counts.all[[puck]][, colSums(counts.all[[puck]]) >= umi.thresh]
}

counts.all <- t(merge.sparse(lapply(counts.all, t)))

# Run Giotto ####
# create instructions
ginstuc <- createGiottoInstructions(python_path = python_path,
                                    show_plot = F, save_plot = F)

ssg <- createGiottoObject(raw_exprs = counts.all, instructions = ginstuc) # slide-seq Giotto


# Filtering and normalization
ssg <- filterGiotto(ssg, expression_values = 'raw', expression_threshold = 1,
                    gene_det_in_min_cells = 150, min_det_genes_per_cell = 0)


# Loading the reference ####
load(seurat.path)
KBR <- subset(KBR, cells = colnames(KBR)[!is.na(KBR@meta.data$subclass.l1)])
KBR <- SetIdent(KBR, value = 'subclass.l1')
KBR@meta.data[['nUMI']] <- KBR@meta.data$nCount_RNA


# removing some unrelated l1 subclasses ####
subc2remove <- c('PapE', "NEU")
KBR.f <- subset(KBR, cells = colnames(KBR)[!KBR@meta.data$subclass.l1 %in% subc2remove])
rm(KBR)

# Removing some low count or injured level2 subclasses ####
subc2remove.l2 <- c('dATL', 'dDTL3', 'dPT', 'dFIB', 'dC-TAL', 'dEC-PTC', 'dVSMC/P',
                    'dC-IC-A', 'dCNT', 'dEC', 'dPOD', 'dDCT',
                    'dM-FIB', 'dIMCD', 'dM-PC',
                    'dPT/DTL', 'dM-TAL', 'dVSMC', 'aTAL2', "B", "N")

KBR.f <- subset(KBR.f, cells = colnames(KBR.f)[!KBR.f@meta.data$subclass.l2 %in% subc2remove.l2])

KBR.f@meta.data[['labels.rctd.l1']] <- KBR.f@meta.data[['subclass.l1']]
KBR.f@meta.data[['labels.rctd.l1']][KBR.f@meta.data[['labels.rctd.l1']] %in% c("DTL", "ATL")] <- "TL" # joining the thin limbs at level1
KBR.f@meta.data[['labels.rctd.l1']] <- gsub("/", '_', KBR.f@meta.data[['labels.rctd.l1']]) # RCTD doesn't like /

KBR.f@meta.data[['labels.rctd.l2']] <- KBR.f@meta.data[['subclass.l2']]
KBR.f@meta.data[['labels.rctd.l2']] <- gsub("/", '_', KBR.f@meta.data[['labels.rctd.l2']]) # RCTD doesn't like /

KBR.f <- SetIdent(KBR.f, value = 'labels.rctd.l2')
KBR.f <- subset(KBR.f, downsample = 3000)

# running l1.RCTD ####
print('Running RCTD on subclass.l1')
celltypes.rctd <- KBR.f@meta.data[['labels.rctd.l1']]
celltypes.rctd <- factor(celltypes.rctd, levels = sort(unique(celltypes.rctd)))
names(celltypes.rctd) <- rownames(KBR.f@meta.data)

numi.rctd <- setNames(KBR.f@meta.data[['nUMI']], rownames(KBR.f@meta.data))
sn.ref <- Reference(counts = KBR.f@assays$RNA@counts, cell_types = celltypes.rctd, nUMI = numi.rctd)
spat.rna <- SpatialRNA(counts = ssg@raw_exprs, use_fake_coords = T)
rctd.l1 <- create.RCTD(spatialRNA = spat.rna, reference = sn.ref, max_cores = 6,
                        gene_cutoff = 0.0003, gene_cutoff_reg = 0.00035, fc_cutoff = 0.45, fc_cutoff_reg = 0.6)
if(!"REN" %in% rctd.l1@internal_vars$gene_list_reg)
  rctd.l1@internal_vars$gene_list_reg <- c(rctd.l1@internal_vars$gene_list_reg, "REN") # adding REN

if(!"REN" %in% rctd.l1@internal_vars$gene_list_bulk)
  rctd.l1@internal_vars$gene_list_bulk <- c(rctd.l1@internal_vars$gene_list_bulk, "REN") # adding REN

rctd.l1@spatialRNA@counts <- rbind(rctd.l1@spatialRNA@counts, 
                                   REN=ssg@raw_exprs['REN', colnames(rctd.l1@spatialRNA@counts)])

rctd.l1 <- fitBulk(rctd.l1)
rctd.l1 <- choose_sigma_c(rctd.l1)
rctd.l1 <- fitPixels_full(rctd.l1)
colnames(rctd.l1@results$weights) <- gsub("_", "/", colnames(rctd.l1@results$weights)) # revert _ to / in celltype names
rctd.l1.list <- list(config = rctd.l1@config, cell_type_info = rctd.l1@cell_type_info, results = rctd.l1@results,
                      ref.cellIDs = colnames(rctd.l1@reference@counts), ref.genes = rownames(rctd.l1@reference@counts),
                      spRNA.cellIDs = colnames(rctd.l1@spatialRNA@counts), spRNA.genes = rownames(rctd.l1@spatialRNA@counts),
                      internal_vars = rctd.l1@internal_vars)
saveRDS(rctd.l1.list, file = l1file)

# running l2.RCTD ####

# # forming a l1 to l2 transition table
l1.l2 <- as.data.frame(table(gsub('_', '/', KBR.f$labels.rctd.l1), gsub("_", "/", KBR.f$labels.rctd.l2)))
l1.l2 <- l1.l2[l1.l2$Freq > 0, ]
l1.l2 <- l1.l2[order(l1.l2$Var1), ]
colnames(l1.l2) <- c('name.l1', 'name.l2', 'Freq.KBR')
write.table(l1.l2, paste0(rctd.dir, '/', prefix.f, '_l1-l2_subclasses.tsv'), sep = '\t', quote = F, row.names = F)

# running l2.RCTD ####
KBR.f2 <- KBR.f
rm(KBR.f)

# Load the RCTD results of subclass l1
l1.weight.thresh <- 40 # the minimum weight to further decompose in subclass l2
rctd.list <- readRDS(l1file)
ssg <- subsetGiotto(ssg, cell_ids = rctd.list$spRNA.cellIDs)
ssg <- addRCTDMaxWeight(ssg, rctd.list, meta.prefix = 'rctd.l1.')

# merging some subclass.l2 populations
subc2merge <- c('')
subc2merge.dest <- c('')
KBR.f2[['subclass.l2.mod']] <- KBR.f2$labels.rctd.l2
for(i in seq_along(subc2merge)){
  KBR.f2[['subclass.l2.mod']][KBR.f2[['subclass.l2.mod']] == subc2merge[i]] <- subc2merge.dest[i]
}

KBR.f2 <- SetIdent(KBR.f2, value = 'subclass.l2.mod')

# Running RCTD on the subtypes of each l1 subclass
l1.classes <- unique(KBR.f2$labels.rctd.l1)
rctd.list.l2 <- list()
for (l1cl in l1.classes){
  print(paste0('RCTD on ', l1cl))
  print(paste0(Sys.time(), ': ', 'RCTD on ', l1cl))

  # Keeping on l1 subclasses + l2 subclasses of l1cl
  combinedClasses <- KBR.f2$labels.rctd.l1
  combinedClasses[combinedClasses == l1cl] <- KBR.f2$subclass.l2.mod[combinedClasses == l1cl]
  combinedClasses <- factor(combinedClasses, levels = sort(unique(combinedClasses)))

  celltypes.rctd.l2 <- combinedClasses
  names(celltypes.rctd.l2) <- rownames(KBR.f2@meta.data)

  print(paste0("L2 subclasses include: "))
  print(unique(KBR.f2$subclass.l2.mod[KBR.f2$labels.rctd.l1 == l1cl]))

  # beads classified as l1cl on RCTD.l1 run
  ssg.cl <- subsetGiotto(ssg, cell_ids = ssg@cell_ID[ssg@cell_metadata$rctd.l1.maxCelltype == gsub("_", "/", l1cl)])
  ssg.cl <- subsetGiotto(ssg.cl, cell_ids = ssg.cl@cell_ID[ssg.cl@cell_metadata$rctd.l1.maxWeight >= l1.weight.thresh])

  # running RCTD
  numi.rctd <- setNames(KBR.f2@meta.data[['nUMI']], rownames(KBR.f2@meta.data))
  sn.ref <- Reference(counts = KBR.f2@assays$RNA@counts, cell_types = celltypes.rctd.l2, nUMI = numi.rctd)
  spat.rna <- SpatialRNA(counts = ssg@raw_exprs, use_fake_coords = T)
  rctd.cl <- create.RCTD(spatialRNA = spat.rna, reference = sn.ref, max_cores = 6,
                          gene_cutoff = 0.0003, gene_cutoff_reg = 0.00035, fc_cutoff = 0.45, fc_cutoff_reg = 0.6)
  rctd.cl <- fitBulk(rctd.cl)
  rctd.cl <- choose_sigma_c(rctd.cl)
  rctd.cl@spatialRNA@counts <- rctd.cl@spatialRNA@counts[, ssg.cl@cell_ID]
  rctd.cl@spatialRNA@nUMI <- rctd.cl@spatialRNA@nUMI[ssg.cl@cell_ID]
  rctd.cl@spatialRNA@coords <- rctd.cl@spatialRNA@coords[ssg.cl@cell_ID, ]
  rctd.cl <- fitPixels_full(rctd.cl)
  colnames(rctd.cl@results$weights) <- gsub("_", "/", colnames(rctd.cl@results$weights)) # revert _ to / in celltype names
  rctd.cl.list <- list(config = rctd.cl@config, cell_type_info = rctd.cl@cell_type_info, results = rctd.cl@results,
                       ref.cellIDs = colnames(rctd.cl@reference@counts), ref.genes = rownames(rctd.cl@reference@counts),
                       spRNA.cellIDs = colnames(rctd.cl@spatialRNA@counts), spRNA.genes = rownames(rctd.cl@spatialRNA@counts),
                       internal_vars = rctd.cl@internal_vars)

  rctd.list.l2[[gsub("_", "/", l1cl)]] <- rctd.cl.list
  saveRDS(rctd.list.l2, file = l2file)
  print(paste0(Sys.time(), ": names in rctd.list.l2:"))
  print(names(rctd.list.l2))
}
sink()


# Combining l1 and l2 weights into a single matrix####
rctd.l1 <- readRDS(l1file)
rctd.l2 <- readRDS(l2file)
l1.l2 <- read.table(paste0(rctd.dir, '/', prefix.f, '_l1-l2_subclasses.tsv'), sep = '\t', header = T)
weights <- list(rctd.l1$results$weights)
weights[[1]] <- 100 * sweep(weights[[1]], 1, rowSums(weights[[1]]), "/")
colnames(weights[[1]]) <- paste0('l1.', colnames(weights[[1]]))

for(l1cl in unique(l1.l2$name.l1)){
  l2cls <- unique(l1.l2$name.l2[l1.l2$name.l1 == l1cl])
  rctd.clsl2 <- rctd.l2[[l1cl]]$results$weights
  rctd.clsl2 <- 100 * sweep(rctd.clsl2, 1, rowSums(rctd.clsl2), "/")
  colnames(rctd.clsl2) <- paste0('l2.', colnames(rctd.clsl2))
  weights <- append(weights, rctd.clsl2)
}

all.colnames <- unique(unlist(lapply(weights, colnames)))
weights.all <- data.frame(matrix(nrow = dim(weights[[1]])[1],
                                 ncol = length(all.colnames)))
rownames(weights.all) <- rownames(weights[[1]])
colnames(weights.all) <- all.colnames
for(mat in weights){
  weights.all[rownames(mat), colnames(mat)] <- mat
}
weights.all[is.na(weights.all)] <- 0
write.table(weights.all,
            file = paste0(rctd.dir, '/', prefix.f, '_combinedWeights_norm.tsv'),
            sep = '\t', quote = F)


# Save all pucks as Giotto objects, including the boundaries if exists
weights.all.norm <- read.table(paste0(rctd.dir, '/', prefix.f, '_combinedWeights_norm.tsv'), sep = '\t', check.names = F)

foreach(puck=names(puckdirs)) %do%{
  ssg <- readPuck_j(puck, datadirs,
                    weights.all.norm, instrs = ginstuc)
  ssg <- addCellMetadata(ssg, new_metadata = colSums(ssg@raw_exprs), vector_name = 'nUMI')
  ssg <- addCellMetadata(ssg, IsBetweenLines(ssg, ex.lines[[puck]]), vector_name = 'inBoundary')
  if(sum(is.na(ssg@cell_metadata$inBoundary)) > 0 ){
    print("problem with puck", puck)
    return()
  }
  saveRDS(ssg, file = paste0(gioObjFolder, "/", puck, ".rds"))
  puck
}


# Plot the l1 weights ####
plot_f <- paste0(plot.dir, "/l1weights")
if(!dir.exists(plot_f))
  dir.create(plot_f)
weights.all <- read.table(paste0(rctd.dir, '/', prefix.f, '_combinedWeights_norm.tsv'), check.names = F)

l1.classes <- colnames(weights.all)[startsWith(colnames(weights.all), 'l1.')]
weight_cutoff <- 30
# for(puck in names(puckdirs)){
foreach(puck=names(puckdirs)) %do%{
  ssg.puck <- readRDS(paste0(gioObjFolder, "/", puck, ".rds"))
  pdf(paste0(plot_f, '/', puck, '-l1weights.pdf'))
  for(l1cl in l1.classes){
    metaCol <- l1cl # paste0("l1.", l1cl)
    sel.cells <- ssg.puck@cell_ID[(ssg.puck@cell_metadata[[metaCol]] >= weight_cutoff) & ssg.puck@cell_metadata[['inBoundary']]]
    plt <- mySpatPlot1(ssg.puck,
                       select_cells = sel.cells,
                       cell_color = metaCol, color_as_factor = F,
                       cell_color_gradient = c('green3', 'orange1', 'orangered4'),
                       gradient_limits = c(20, 80), gradient_midpoint = 50,
                       point_size = 0.2, #other_cells_alpha = 0.1,
                       other_cell_color = 'gray94', other_point_size = 0.1,
                       show_legend = T, ctheme = theme_noAxis(T)) #+ coord_fixed()
    plt <- plt + ggtitle(l1cl) + theme(plot.title = element_text(face = 'bold', vjust = -3, size = 13))
    print(plt)
  }
  dev.off()
}

# Plot the l2 weights ####
plot_f <- paste0(plot.dir, "/l2weights")
if(!dir.exists(plot_f))
  dir.create(plot_f)

l1.l2 <- read.table(paste0(rctd.dir, '/', prefix.f, '_l1-l2_subclasses.tsv'), sep = '\t', header = T)
weights.all <- read.table(paste0(rctd.dir, '/', prefix.f, '_combinedWeights_norm.tsv'), check.names = F)
weights.all[is.na(weights.all)] <- 0

l2.classes <- colnames(weights.all)[startsWith(colnames(weights.all), 'l2.')]
weight_cutoff <- 30
# for(puck in names(puckdirs)){
foreach(puck=names(puckdirs)) %do%{
  ssg.puck <- readRDS(paste0(gioObjFolder, "/", puck, ".rds"))
  pdf(paste0(plot_f, '/', puck, '-l2weights.pdf'))
  for(l2cl in l2.classes){
    metaCol <- l2cl # paste0("l2.", l2cl)
    sel.cells <- ssg.puck@cell_ID[(ssg.puck@cell_metadata[[metaCol]] >= weight_cutoff) & ssg.puck@cell_metadata[['inBoundary']]]
    plt <- mySpatPlot1(ssg.puck,
                       select_cells = sel.cells,
                       cell_color = metaCol, color_as_factor = F,
                       cell_color_gradient = c('green3', 'orange1', 'orangered4'),
                       gradient_limits = c(20, 80), gradient_midpoint = 50,
                       point_size = 0.2, #other_cells_alpha = 0.1,
                       other_cell_color = 'gray94', other_point_size = 0.1,
                       show_legend = T, ctheme = theme_noAxis(T)) #+ coord_fixed()
    plt <- plt + ggtitle(l2cl) + theme(plot.title = element_text(face = 'bold', vjust = -3, size = 13))
    print(plt)
  }
  dev.off()
}


# Plot l1 beads assignment ####
# Set up colors ####
l1.l2 <- read.table(paste0(rctd.dir, '/', prefix.f, '_l1-l2_subclasses.tsv'), sep = '\t', header = T)
weights.all <- read.table(paste0(rctd.dir, '/', prefix.f, '_combinedWeights_norm.tsv'), check.names = F)
l1.colors.df <- colorTable[, c('subclass.l1_label', 'subclass.l1_color')]
l1.colors.df <- l1.colors.df[!duplicated(l1.colors.df), ]
l1.colors.df <- l1.colors.df[l1.colors.df$subclass.l1_label %in% substr(colnames(weights.all), 4, 100)]
l1.colors <- l1.colors.df$subclass.l1_color
names(l1.colors) <- paste0('l1.', l1.colors.df$subclass.l1_label)

l2.colors.df <- colorTable[, c('subclass.l2_label', 'subclass.l2_color')]
l2.colors.df <- l2.colors.df[!duplicated(l2.colors.df), ]
l2.colors.df <- l2.colors.df[l2.colors.df$subclass.l2_label %in% substr(colnames(weights.all), 4, 100)]
l2.colors <- l2.colors.df$subclass.l2_color
names(l2.colors) <- paste0('l2.', l2.colors.df$subclass.l2_label)

pdf(paste0(plot.dir, "/", prefix.f, "_l1Assignments.pdf"), width = 11, height = 11)
for(puck in names(puckdirs)){
  ssg <- readRDS(paste0(gioObjFolder, "/", puck, ".rds"))
  sel.cells <- ssg@cell_ID[(ssg@cell_metadata$has.l1.rctd) &
                             (ssg@cell_metadata$maxWeight.l1 >= 40) &
                             ssg@cell_metadata[['inBoundary']]]
  plt <- mySpatPlot1(ssg, select_cells = sel.cells,
                     cell_color = 'maxCelltype.l1', show_legend = T,
                     cell_color_code = l1.colors, other_cell_color = 'gray95', other_point_size=0.3) +
    theme_noAxis(T)
  plt <- plt + ggtitle(paste0(puck, ' - ', 'maxCelltype.l1')) +
    theme(plot.title = element_text(face = 'bold', vjust = -3, size = 13))
  print(plt)
}
dev.off()

# Plot l2 beads assignment ####
pdf(paste0(plot.dir, "/", prefix.f, "_l2Assignments.pdf"), width = 11, height = 11)
for(puck in names(puckdirs)){
  ssg <- readRDS(paste0(gioObjFolder, "/", puck, ".rds"))
  sel.cells <- ssg@cell_ID[(ssg@cell_metadata$has.l2.rctd) &
                             (ssg@cell_metadata$maxWeight.l2 >= 30) &
                             (substr(ssg@cell_metadata$maxCelltype.l2, 4, 100) %in% l1.l2$name.l2) &
                             ssg@cell_metadata[['inBoundary']]]
  plt <- mySpatPlot1(ssg, select_cells = sel.cells,
                     cell_color = 'maxCelltype.l2', show_legend = T,
                     cell_color_code = l2.colors, other_cell_color = 'gray95', other_point_size=0.3) +
    theme_noAxis(T)
  plt <- plt + ggtitle(paste0(puck, ' - ', 'maxCelltype.l2')) +
    theme(plot.title = element_text(face = 'bold', vjust = -3, size = 13))
  print(plt)
}
dev.off()