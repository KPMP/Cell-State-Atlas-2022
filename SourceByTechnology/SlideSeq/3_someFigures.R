library(ggplot2)
library(Matrix)
library(Giotto)
library(reshape2)
library(patchwork)
source("./utils.R")

python_path <- "/media/Home_Raid1_Voyager/kian/anaconda3/envs/giotto_200703/bin/python"
out.dir <- "./outputs"
slideseq.dir <- "../ForHubmapCollaborators"
rctd.dir <- paste0(out.dir, "/", "deconv_files")
plot.dir <- paste0(out.dir, "/", "Plots")

pucks <- list.dirs(slideseq.dir, full.names = FALSE)
pucks <- pucks[nchar(pucks) > 0]

good.pucks <- c("Puck_200903_26", "Puck_200903_02", "Puck_200903_01", "Puck_200903_03", 
           "Puck_200903_07", "Puck_200903_05", "Puck_200906_18", "Puck_200906_10",
           "Puck_200903_06", "Puck_200903_25")
weights.all <- read.table(paste0(rctd.dir, "/", "combinedWeights.tsv"), 
                          sep = '\t', header = T, check.names = F)

instrs = createGiottoInstructions(save_dir = plot.dir,
                                  save_plot = TRUE,
                                  show_plot = FALSE,
                                  python_path = python_path, plot_format = 'pdf')

l1.l2 <- read.table(paste0(rctd.dir, "/", "rctd_subclasses.tsv"), sep = '\t', header = T)

# Setup colors ####
colorTable <- data.table::fread("Cluster_Color_Table.tsv", select = 1:20)
l1.colors.df <- colorTable[, c('subclass.l1_label', 'subclass.l1_color')]
l1.colors.df <- l1.colors.df[!duplicated(l1.colors.df), ]
l1.colors <- l1.colors.df$subclass.l1_color
names(l1.colors) <- paste0('l1.', l1.colors.df$subclass.l1_label)

l2.colors.df <- colorTable[, c('subclass.l2_label', 'subclass.l2_color')]
l2.colors.df <- l2.colors.df[!duplicated(l2.colors.df), ]
l2.colors <- l2.colors.df$subclass.l2_color
names(l2.colors) <- paste0('l2.', l2.colors.df$subclass.l2_label)


# RCTD Cell type distribution ####
l1Cols <- colnames(weights.all)[grepl("^l1.", colnames(weights.all))]
l2Cols <- colnames(weights.all)[grepl("^l2.", colnames(weights.all))]

l1classes <- l1.colors.df$subclass.l1_label
l2classes <- l2.colors.df$subclass.l2_label
l2classes <- l2classes[l2classes %in% substr(colnames(weights.all), 4, 100)]

weights.selpucks <- weights.all[substr(rownames(weights.all), 1, nchar(good.pucks[1])) %in% good.pucks, ]

# plotting l2 freqs ####
weights.selpucks[is.na(weights.selpucks)] <- 0
l2.maxWeight <- apply(weights.selpucks[, l2Cols], 1, max)
l2.maxWeight.ct <- l2Cols[apply(weights.selpucks[, l2Cols], 1, which.max)]

l1_to_l2_wrong <- !substr(l2.maxWeight.ct, 4, 100) %in% l1.l2$name.l2 # those assignments that don't match between l1 and l2 types
l2.maxWeight <- l2.maxWeight[!l1_to_l2_wrong]
l2.maxWeight.ct <- l2.maxWeight.ct[!l1_to_l2_wrong]


l2.freqs <- as.data.frame(table(l2.maxWeight.ct[l2.maxWeight >= 50]), stringsAsFactors = F)
colnames(l2.freqs) <- c('Type', 'Freq')
l2.freqs$Freq <- l2.freqs$Freq / sum(l2.freqs$Freq)
l2.freqs$color <- l2.colors[l2.freqs$Type]
l2.freqs$Type <- substr(l2.freqs$Type, 4, stop = 100)
l2.freqs$Type <- factor(l2.freqs$Type, levels = l2classes)
pdf(paste0(plot.dir, "/", "ExtFig6a_l2freqs.pdf"), height = 7, width = 11)
ggplot(l2.freqs) + geom_bar(aes(x = Type, y = Freq, fill = Type), stat = 'identity') + 
  scale_fill_manual(values = setNames(l2.freqs[,'color'], l2.freqs$Type)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()




# UMI counts and weights per cell type subclass l2 ####
# finding the UMI counts for all beads of all pucks
umi.counts <- list()
have.l2s <- list()
for(puck in pucks){
  ssg <- readPuck(puck, paste(slideseq.dir, puck, sep = '/'), weights.all, instrs)
  umi.counts[[puck]] <- colSums(ssg@raw_exprs)
  have.l2s[[puck]] <- ssg@cell_metadata$has.l2.rctd
  names(have.l2s[[puck]]) <- ssg@cell_metadata$cell_ID
}
l2Cols <- colnames(weights.all)[grepl("^l2.", colnames(weights.all))]

l2classes <- l2.colors.df$subclass.l2_label
l2classes <- l2classes[l2classes %in% substr(colnames(weights.all), 4, 100)]

# selecting beads that have l2 deconvolution
names(have.l2s) <- NULL
have.l2s <- unlist(have.l2s)
have.l2s <- names(have.l2s)[have.l2s]

weights.selpucks <- weights.all[substr(rownames(weights.all), 1, nchar(good.pucks[1])) %in% good.pucks, ]
weights.selpucks <- weights.selpucks[rownames(weights.selpucks) %in% have.l2s, ]
weights.selpucks[is.na(weights.selpucks)] <- 0

# plotting l2 freqs
l2.maxWeight <- apply(weights.selpucks[, l2Cols], 1, max)
l2.maxWeight.ct <- l2Cols[apply(weights.selpucks[, l2Cols], 1, which.max)]
names(l2.maxWeight.ct) <- rownames(weights.selpucks)

umi.counts.gp <- c() # umi counts of the good pucks
for(gp in good.pucks)
  umi.counts.gp <- append(umi.counts.gp, umi.counts[[gp]])

allumis <- names(umi.counts.gp)
beadsWithWeights <- names(l2.maxWeight)
umisWithWeights <- umi.counts.gp[beadsWithWeights] # UMI counts for those beads that were not filtered

l2.umi.ct <- data.frame(name=beadsWithWeights, UMI=umisWithWeights, maxWeight=l2.maxWeight, Celltype=l2.maxWeight.ct)
l2.umi.ct <- l2.umi.ct[l2.umi.ct$maxWeight >= 50, ]
l2.umi.ct$Celltype <- substr(l2.umi.ct$Celltype, 4, 100)
l2.umi.ct$Celltype <- factor(l2.umi.ct$Celltype, levels = l2classes)

# removing those beads that switched classes in the l2 deconvolution (should be low in number)
l2.umi.ct <- l2.umi.ct[l2.umi.ct$Celltype %in% l1.l2$name.l2, ]

pdf(paste0(plot.dir, "/", "ExtFig6a_UMI_count_per_l2.pdf"), height=6, width = 14)
ggplot(l2.umi.ct, aes(x=Celltype, y = log10(UMI), fill = Celltype)) + 
  geom_violin(draw_quantiles = c(0.5), scale = 'width') + theme_light() + 
  scale_fill_manual(values = setNames(l2.colors, substr(names(l2.colors), 4, 100))) + 
  ggtitle("UMI count per subclass l2 celltypes") + theme(plot.title = element_text(face='bold', hjust=0.5)) + 
  theme(axis.text.x = element_text(angle=90, hjust = 1))
dev.off()


# Violin plots for the weights of called beads
pdf(paste0(plot.dir, "/", "ExtFig6a_weight_per_l2.pdf"), height=6, width = 14)
ggplot(l2.umi.ct, aes(x=Celltype, y = maxWeight, fill = Celltype)) + 
  geom_violin(draw_quantiles = c(0.5), scale = 'width') + theme_light() + 
  scale_fill_manual(values = setNames(l2.colors, substr(names(l2.colors), 4, 100))) + 
  ggtitle("RCTD weight per subclass l2 celltypes") + theme(plot.title = element_text(face='bold', hjust=0.5)) + 
  theme(axis.text.x = element_text(angle=90, hjust = 1))
dev.off()




# Marker gene dotplot for good pucks ####
require(Seurat)
weights.all <- read.table(paste0(rctd.dir, "/", "combinedWeights.tsv"), 
                          sep = '\t', header = T, check.names = F)

l2Cols <- colnames(weights.all)[grepl("^l2.", colnames(weights.all))]

l2classes <- l2.colors.df$subclass.l2_label
l2classes <- l2classes[l2classes %in% substr(colnames(weights.all), 4, 100)]

ssg.list <- list()
for(puck in good.pucks){
  ssg <- readPuck(puck, paste(slideseq.dir, puck, sep = '/'), weights.all, instrs)
  ssg <- subsetGiotto(ssg, cell_ids = ssg@cell_ID[ssg@cell_metadata$maxWeight.l2 >= 50])
  ssg <- subsetGiotto(ssg, cell_ids = ssg@cell_ID[substr(ssg@cell_metadata$maxCelltype.l2, 4, 100) %in% l1.l2$name.l2])
  ssg.list[[puck]] <- ssg
}

counts.all <- t(merge.sparse(lapply(lapply(ssg.list, function(ss) ss@raw_exprs), t)))
meta.all <- do.call(rbind, lapply(ssg.list, function(ss) ss@cell_metadata))
ssg.all <- createGiottoObject(raw_exprs = counts.all, cell_metadata = meta.all)

genes2plot <- c("NPHS2", "PODXL"  , "CFH"   , "VCAM1" , "PRODH2", "GPX3"  , "SLC34A1" , "ALDOB" , "ACSM2B", "ITGB8"  , "SLC12A1",
                "UMOD" , "EGF"    , "PROM1" , "DCDC2" , "ITGB6" , "NOS1"  , "SLC12A3" , "FXYD2" , "SLC8A1", "HSD11B2", "AQP2"   ,
                "AQP3" , "SLC26A7", "SLC4A9", "EMCN"  , "IGFBP5", "TM4SF1", "DNASE1L3", "PDGFRB", "REN"   , "MYH11"  , "MYL9"   ,
                "RGS5" , "COL1A1" , "ACTA2" , "C7"    , "FGF14" , "IGF1"  , "PTPRC"   , "IL7R"  , "CD96"  , "GNLY"   , "GZMA"   ,
                "STAB1", "CSF2RA" , "FLT3"  , "S100A9")

# changing the order C-TAL and aTAL1
l2classes2 <- l2classes
l2classes2[l2classes == 'C-TAL'] <- 'aTAL1'
l2classes2[l2classes == 'aTAL1'] <- 'C-TAL'

rownames(ssg.all@cell_metadata) <- ssg.all@cell_metadata$cell_ID
sss <- CreateSeuratObject(counts = ssg.all@raw_exprs, meta.data = ssg.all@cell_metadata)
sss$maxCelltype.l2 <- factor(substr(sss$maxCelltype.l2, 4, 100), levels = l2classes2)
sss <- SetIdent(sss, value = sss$maxCelltype.l2)

pdf(paste0(plot.dir, "/", "ExtFig6b_markers.pdf"), height=8, width = 12)
DotPlot(sss, features = (genes2plot), cols = c('gray95', 'orange'), ) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5, size = 9.4))
dev.off()



# plotting beads with subclass l2 annotations ####
weights.all <- read.table(paste0(rctd.dir, "/", "combinedWeights.tsv"), 
                          sep = '\t', header = T, check.names = F)
l1.l2 <- read.table(paste0(rctd.dir, "/", "rctd_subclasses.tsv"), sep = '\t', header = T)

pdf(plot.dir, "/", "ExtFig6c_l2assignment.pdf", width = 11, height = 11)
for(puck in pucks){
  ssg.puck <- readPuck(puck, paste(slideseq.dir, puck, sep = '/'), weights.all, instrs)
  
  sel.cells <- ssg.puck@cell_ID[(ssg.puck@cell_metadata$has.l2.rctd) & 
                                  (ssg.puck@cell_metadata$maxWeight.l2 >= 50) & 
                                  (substr(ssg.puck@cell_metadata$maxCelltype.l2, 4, 100) %in% l1.l2$name.l2)]
  plt <- mySpatPlot1(ssg.puck, select_cells = sel.cells, 
                     cell_color = 'maxCelltype.l2', show_legend = T, 
                     cell_color_code = l2.colors, other_cell_color = 'gray95') +
    theme_noAxis(T) 
  plt <- plt + ggtitle(paste0(puck, ' - ', 'maxCelltype.l2')) + 
    theme(plot.title = element_text(face = 'bold', vjust = -3, size = 13))
  print(plt)
}
dev.off()


# # Plotting specific genes ####
# puck <- 'Puck_200903_06'
# ssg.puck <- readPuck(puck, paste(slideseq.dir, puck, sep = '/'), weights.all, instrs)
# plt1 <- plotGenesRGB3(ssg.puck, genes = c("REN", "NOS1", "NPHS2"), max.expr = c(7, 13, 7), legend = 'right',
#                       color.hues = c(0, 0.33, 0.66))
# plt2 <- plotGenesRGB3(ssg.puck, genes = c("MYH11", "SLC5A12", "EMCN"), max.expr = c(10, 5, 6), legend = 'right',
#                       color.hues = c(0, 0.33, 0.66))
# pdf(paste0(plot.dir, "/", "Fig2f_g.pdf"), width = 12, height = 11)
# print(plt1)
# print(plt2)
# dev.off()


# puck <- 'Puck_200903_01'
# ssg.puck <- readPuck(puck, paste(slideseq.dir, puck, sep = '/'), weights.all, instrs)
# plt3 <- plotGenesRGB3(ssg.puck, genes = c("ITGB6", "EGF", "SLC12A1"), max.expr = c(6, 8, 4), legend = 'bottom',
#                       color.hues = c(0, 0.33, 0.66))
# plt4 <- plotMetadataRGB(ssg.puck, metaCols = c("l2.MYOF", "l2.aFIB", "l2.T"), min.expr = c(20, 20, 20),
#                         max.expr = c(80, 55, 70), legend = 'bottom', color.hues = c(0, 0.11, 0.375))
# pdf(paste0(plot.dir, "/", "Fig3c_e.pdf"), width = 12, height = 11)
# print(plt3)
# print(plt4)
# dev.off()