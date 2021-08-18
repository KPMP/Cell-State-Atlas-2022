# How to point at the data
# How to point at the 10x Seurat object

library(ggplot2)
library(Matrix)
library(Giotto)
library(reshape2)
library(patchwork)
library(RCTD)
library(Seurat)
source("./utils.R")

createSpatRNA <- function(gObj, fakeCoord = F){
  coords <- as.data.frame(gObj@spatial_locs)
  colnames(coords) <- c('x', 'y', 'z', 'cells')
  rownames(coords) <- coords[, 'cells']
  coords <- coords[, c('x', 'y')]
  
  counts <- as(gObj@raw_exprs, 'dgCMatrix')
  spat.RNA <- RCTD:::SpatialRNA(coords = coords, counts = counts, nUMI = colSums(counts))
  return(spat.RNA)
}

slideseq.dir <- "../ForHubmapCollaborators"
python_path <- "/media/Home_Raid1_Voyager/kian/anaconda3/envs/giotto_200703/bin/python"
seurat.path <- "../SeuratObjects/Kidney_KPMP-Biopsy_Ref_10X-R_092020_Seurat_F.rda"
  
out.dir <- "./outputs"
if (!dir.exists(out.dir))
  dir.create(out.dir)

rctd.dir <- paste0(out.dir, "/", "deconv_files")
if (!dir.exists(rctd.dir))
  dir.create(rctd.dir)

plot.dir <- paste0(out.dir, "/", "Plots")
if (!dir.exists(plot.dir))
  dir.create(plot.dir)

## Log file
sink(paste0(out.dir, "/", "rctd_run.log"))

# Reading all the pucks at once, making one Giotto object to contain them all, discarding the spatial locations
# Reading counts ####
umi.thresh <- 80
pucks <- list.dirs(slideseq.dir, full.names = FALSE)
pucks <- pucks[nchar(pucks) > 0]
counts.all <- list()
for(puck in pucks){
  print(paste('Reading', puck))
  puck.addr <- paste(slideseq.dir, puck, sep = '/')
  counts.all[[puck]] <- getPuckCounts(puck, puck.addr)
  counts.all[[puck]] <- counts.all[[puck]][, colSums(counts.all[[puck]]) >= umi.thresh]
}

counts.all <- t(merge.sparse(lapply(counts.all, t)))

# Run Giotto ####
# create instructions
instrs = createGiottoInstructions(save_dir = plot.dir,
                                  save_plot = FALSE,
                                  show_plot = FALSE,
                                  python_path = python_path, plot_format = 'pdf')

ssg <- createGiottoObject(raw_exprs = counts.all, instructions = instrs) # slide-seq Giotto

# Filtering and normalization 
ssg <- filterGiotto(ssg, expression_values = 'raw', expression_threshold = 1, 
                    gene_det_in_min_cells = 150, min_det_genes_per_cell = 0)


# Loading the reference ####
load(seurat.path)
KBR <- subset(KBR, cells = colnames(KBR)[!is.na(KBR@meta.data$subclass.l1)])
KBR <- SetIdent(KBR, value = 'subclass.l1')
KBR@meta.data[['nUMI']] <- KBR@meta.data$nCount_RNA


# removing some medullary subclasses ####
subc2remove <- c('ATL', 'DTL', 'NEU', 'PapE')
KBR.f <- subset(KBR, cells = colnames(KBR)[!KBR@meta.data$subclass.l1 %in% subc2remove])
KBR.f@meta.data[['liger_ident_coarse']] <- KBR.f@meta.data[['subclass.l1']]
KBR.f@meta.data[['liger_ident_coarse']] <- factor(KBR.f@meta.data[['liger_ident_coarse']], 
                                                  levels = sort(unique(KBR.f@meta.data[['liger_ident_coarse']])))

# Removing some low count or medullary level2 subclasses ####
subc2remove.l2 <- c('dPT', 'dFIB', 'dC-TAL', 'dEC-PTC', 'dVSMC/P', 
                    'dC-IC-A', 'dCNT', 'dEC', 'dPOD', 'dDCT',
                    'EC-AVR', 'EC-DVR', 'EC-LYM', 
                    'dM-FIB', 'M-FIB', 'M-IC-A', 
                    'tPC-IC', 'MAC-M2', 'MDC', 'ncMON', 'pDC',
                    'MAST', 'PL', 'dIMCD', 'dM-PC', 'IMCD', 'M-PC',
                    'dPT/DTL', 'dM-TAL', 'M-TAL', 'dVSMC', 'aTAL2')


KBR.f <- subset(KBR.f, cells = colnames(KBR.f)[!KBR.f@meta.data$subclass.l2 %in% subc2remove.l2])
rm(KBR)


# running l1.RCTD ####
print('Running RCTD on subclass.l1')
spat.rna <- createSpatRNA(ssg, fakeCoord = T)
rctd.all <- create.RCTD(spatialRNA = spat.rna, reference = KBR.f, max_cores = 8,
                        gene_cutoff = 0.0001, gene_cutoff_reg = 0.00015, fc_cutoff = 0.4, fc_cutoff_reg = 0.5)
rctd.all <- run.RCTD(rctd.all, doublet_mode = TRUE)

rctd.all.list <- list(config = rctd.all@config, cell_type_info = rctd.all@cell_type_info, results = rctd.all@results,
                      ref.cellIDs = colnames(rctd.all@reference), ref.genes = rownames(rctd.all@reference),
                      spRNA.cellIDs = colnames(rctd.all@spatialRNA@counts), spRNA.genes = rownames(rctd.all@spatialRNA@counts),
                      internal_vars = rctd.all@internal_vars)
saveRDS(rctd.all.list, file = paste0(rctd.dir, "/", "AllPucks_l1results.rds"))


# running l2.RCTD ####
KBR.f2 <- KBR.f
rm(KBR.f)

# Load the RCTD results of subclass l1 
l1.weight.thresh <- 50 # the minimum weight to further decompose in subclass l2
rctd.list <- readRDS(paste0(rctd.dir, "/", "AllPucks_l1results.rds"))
ssg <- subsetGiotto(ssg, cell_ids = rctd.list$spRNA.cellIDs)
ssg <- addRCTDMaxWeight(ssg, rctd.list, meta.prefix = 'rctd.l1.')

# merging some subclass.l2 populations 
subc2merge <- c('')
subc2merge.dest <- c('')
KBR.f2[['subclass.l2.mod']] <- KBR.f2$subclass.l2
for(i in seq_along(subc2merge)){
  KBR.f2[['subclass.l2.mod']][KBR.f2[['subclass.l2.mod']] == subc2merge[i]] <- subc2merge.dest[i]
}

# # Saving the UMAP that includes the l2 subtypes
# KBR.f2 <- SetIdent(KBR.f2, value = 'subclass.l2.mod')
# pdf(paste0(plot.dir, "/", '210428_RCTD_BulkCorrection_UMAP.pdf'), width = 11)
# DimPlot(KBR.f2, label = T)
# dev.off()


# forming a l1 to l2 transition table
l1.l2 <- as.data.frame(table(KBR.f2$subclass.l1, KBR.f2$subclass.l2.mod))
l1.l2 <- l1.l2[l1.l2$Freq > 0, ]
l1.l2 <- l1.l2[order(l1.l2$Var1), ]
colnames(l1.l2) <- c('name.l1', 'name.l2', 'Freq.KBR')
write.table(l1.l2, paste0(rctd.dir, "/", "rctd_subclasses.tsv"), sep = '\t', quote = F, row.names = F)

# Running RCTD on the subtypes of each l1 subclass 
l1.classes <- unique(l1.l2$name.l1)
rctd.list.l2 <- list()
for (l1cl in l1.classes){
  print(paste0('RCTD on ', l1cl))
  print(paste0(Sys.time(), ': ', 'RCTD on ', l1cl))
  
  # Keeping on l1 subclasses + l2 subclasses of l1cl
  combinedClasses <- KBR.f2$subclass.l1
  combinedClasses[combinedClasses == l1cl] <- KBR.f2$subclass.l2.mod[combinedClasses == l1cl]
  combinedClasses <- factor(combinedClasses, levels = sort(unique(combinedClasses)))
  KBR.f2[['liger_ident_coarse']] <- combinedClasses
  
  print(paste0("L2 subclasses include: "))
  print(unique(KBR.f2$subclass.l2.mod[KBR.f2$subclass.l1 == l1cl]))
  
  print(paste0("L2 subclasses include: "))
  print(unique(KBR.f2$subclass.l2.mod[KBR.f2$subclass.l1 == l1cl]))
  
  # beads classified as l1cl on RCTD.l1 run
  ssg.cl <- subsetGiotto(ssg, cell_ids = ssg@cell_ID[ssg@cell_metadata$rctd.l1.maxCelltype == l1cl])
  ssg.cl <- subsetGiotto(ssg.cl, cell_ids = ssg.cl@cell_ID[ssg.cl@cell_metadata$rctd.l1.maxWeight >= l1.weight.thresh])
  
  # running RCTD
  spat.rna <- createSpatRNA(ssg, fakeCoord = T)
  rctd.cl <- create.RCTD(spatialRNA = spat.rna, reference = KBR.f2, max_cores = 8,
                         gene_cutoff = 0.0001, gene_cutoff_reg = 0.00015, 
                         fc_cutoff = 0.4, fc_cutoff_reg = 0.5, test_mode = F)
  rctd.cl <- fitBulk(rctd.cl)
  rctd.cl <- choose_sigma_c(rctd.cl)
  rctd.cl@spatialRNA@counts <- rctd.cl@spatialRNA@counts[, ssg.cl@cell_ID]
  rctd.cl@spatialRNA@nUMI <- rctd.cl@spatialRNA@nUMI[ssg.cl@cell_ID]
  rctd.cl@spatialRNA@coords <- rctd.cl@spatialRNA@coords[ssg.cl@cell_ID, ]
  rctd.cl <- fitPixels(rctd.cl, doublet_mode = TRUE)
  
  rctd.cl.list <- list(config = rctd.cl@config, cell_type_info = rctd.cl@cell_type_info, results = rctd.cl@results,
                       ref.cellIDs = colnames(rctd.cl@reference), ref.genes = rownames(rctd.cl@reference),
                       spRNA.cellIDs = colnames(rctd.cl@spatialRNA@counts), spRNA.genes = rownames(rctd.cl@spatialRNA@counts),
                       internal_vars = rctd.cl@internal_vars)
  rctd.list.l2[[l1cl]] <- rctd.cl.list
  saveRDS(rctd.list.l2, file = paste0(rctd.dir, "/", "AllPucks_l2results.rds"))
  
  print(paste0(Sys.time(), ": names in rctd.list.l2:"))
  print(names(rctd.list.l2))
}
sink()


# combining weights from l1.RCTD and l2.RCTD into a single matrix
rctd.l1 <- readRDS(paste0(rctd.dir, "/", "AllPucks_l1results.rds"))
rctd.l2 <- readRDS(paste0(rctd.dir, "/", "AllPucks_l2results.rds"))
l1.l2 <- read.table(paste0(rctd.dir, "/", "rctd_subclasses.tsv"), sep = '\t', header = T)
weights <- list(rctd.l1$results$weights)
weights[[1]] <- 100 * weights[[1]] / rowSums(weights[[1]])
colnames(weights[[1]]) <- paste0('l1.', colnames(weights[[1]]))

for(l1cl in unique(l1.l2$name.l1)){
  l2cls <- unique(l1.l2$name.l2[l1.l2$name.l1 == l1cl])
  rctd.clsl2 <- rctd.l2[[l1cl]]$results$weights
  rctd.clsl2 <- 100 * rctd.clsl2 / rowSums(rctd.clsl2)
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
write.table(weights.all, file = paste0(rctd.dir, "/", "combinedWeights.tsv"), 
            sep = '\t', quote = F)
