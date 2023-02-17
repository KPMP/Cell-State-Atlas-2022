# Slide-seq - Figure panels (Fig 2,5, ED Fig 4,11) 
library(ggplot2)
library(Matrix)
library(Giotto)
library(reshape2)
library(patchwork)
library(Seurat)
library(doParallel)
library(foreach)
library(Seurat)
library(dplyr)
library(tidyr)
registerDoParallel(5)

source("GiottoUtils.R")
colorTable <- readRDS("Cluster.Color.Table.RDS")
python_path <- "/anaconda3/envs/giotto_200703/bin/python"

#Giotto Objects 
gioObjFolders <- "/GiottoObjects_norm/"
l1.l2 <- read.table("l1-l2_subclasses.tsv", header = T)

# prepping puck names and addressed
pucks <- foreach(gDir = gioObjFolders) %do% {
  pucks <- list.files(gDir)
  pucks <- pucks[startsWith(pucks, "Puck_") & endsWith(pucks, ".rds")]
  pucks <- substr(pucks, 1, nchar(pucks) - 4)
  pucks
}

# setting up the color table
l1.colors.df <- colorTable[, c('subclass.l1_label', 'subclass.l1_color')]
l1.colors.df <- l1.colors.df[!duplicated(l1.colors.df), ]
l1.colors.df <- l1.colors.df[l1.colors.df$subclass.l1_label %in% l1.l2$name.l1]
l1.colors <- l1.colors.df$subclass.l1_color
names(l1.colors) <- paste0('l1.', l1.colors.df$subclass.l1_label)

l2.colors.df <- colorTable[, c('subclass.l2_label', 'subclass.l2_color')]
l2.colors.df <- l2.colors.df[!duplicated(l2.colors.df), ]
l2.colors.df <- l2.colors.df[l2.colors.df$subclass.l2_label %in% l1.l2$name.l2]
l2.colors <- l2.colors.df$subclass.l2_color
names(l2.colors) <- paste0('l2.', l2.colors.df$subclass.l2_label)

# Giotto's parameters
ginstuc <- createGiottoInstructions(python_path = python_path,
                                    show_plot = F, save_plot = F)

#Load Giotto objects for pucks of interest
pucks.use <- c("Puck_200903_01","Puck_200903_06", "Puck_210119_09", "Puck_210119_18")
ssg.list <- foreach(puck=pucks.use, .final=function(x) setNames(x, pucks.use)) %dopar% {
  i <- if(puck %in% pucks[[1]]) 1 else 2
  ssg <- readRDS(paste0(gioObjFolders[[i]], puck, ".rds"))
  ssg <- subsetGiotto(ssg, cell_ids = ssg@cell_ID[ssg@cell_metadata$has.l1.rctd])
  ssg <- subsetGiotto(ssg, cell_ids = ssg@cell_ID[ssg@cell_metadata$inBoundary])
  ssg 
}







###Plot ED Fig 4c
ssg.puck <- ssg.list[[1]]
ssg.puck <- normalizeGiotto(ssg.puck, scalefactor = 10000, log_norm = T, scale_genes = T, scale_cells = F)
puck = "200903_01"

#Remove outlier beads
spatPlot(ssg.puck)
ssg.puck <- subsetGiottoLocs(ssg.puck, x_min=800, x_max=5700, y_min=300, y_max=5100)

#Plot subclass.l2
sel.cells <- ssg.puck@cell_ID[(ssg.puck@cell_metadata$has.l2.rctd) & 
                                (ssg.puck@cell_metadata$maxWeight.l2 >= 30) & 
                                (substr(ssg.puck@cell_metadata$maxCelltype.l2, 4, 100) %in% l1.l2$name.l2)]
plt3 <- mySpatPlot1(ssg.puck, select_cells = sel.cells, 
                    cell_color = 'maxCelltype.l2', show_legend = T, 
                    cell_color_code = l2.colors, other_cell_color = 'gray95') +
  theme_noAxis(T) 
plt3 <- plt3 + ggtitle(paste0(puck, ' - ', 'maxCelltype.l2')) + 
  theme(plot.title = element_text(face = 'bold', vjust = -3, size = 13))
print(plt3)

#Plot Gloms
to.use <- c("l2.POD","l2.PEC","l2.MC","l2.REN","l2.EC-AEA","l2.EC-GC")
ssg.puck.f <- subsetGiotto(ssg.puck, ssg.puck@cell_ID[ssg.puck@cell_metadata$maxCelltype.l2 %in% to.use])

sel.cells <- ssg.puck.f@cell_ID[(ssg.puck.f@cell_metadata$has.l2.rctd) & 
                                  (ssg.puck.f@cell_metadata$maxWeight.l2 >= 30) & 
                                  (substr(ssg.puck.f@cell_metadata$maxCelltype.l2, 4, 100) %in% l1.l2$name.l2)]
plt5 <- mySpatPlot1(ssg.puck.f, select_cells = sel.cells, 
                    cell_color = 'maxCelltype.l2', show_legend = T, 
                    cell_color_code = l2.colors, other_cell_color = 'gray95') +
  theme_noAxis(T) 
plt5 <- plt5 + ggtitle(paste0(puck, ' - ', 'maxCelltype.l2')) + 
  theme(plot.title = element_text(face = 'bold', vjust = -3, size = 13))
print(plt5)

#Plot tubules
to.use <- c("l2.CNT","l2.M-TAL","l2.M-PC","l2.IMCD","l2.DTL3","l2.DTL1","l2.ATL","l2.DTL2",
            "l2.C-PC","l2.aPT","l2.PT-S3","l2.DCT","l2.aTAL1","l2.C-TAL","l2.MD",
            "l2.PT-S1", "l2.IC-B", "l2.PT-S2","l2.PT-S3","l2.C-IC-A","l2.M-IC-A")
ssg.puck.f <- subsetGiotto(ssg.puck, ssg.puck@cell_ID[ssg.puck@cell_metadata$maxCelltype.l2 %in% to.use])

sel.cells <- ssg.puck.f@cell_ID[(ssg.puck.f@cell_metadata$has.l2.rctd) & 
                                  (ssg.puck.f@cell_metadata$maxWeight.l2 >= 30) & 
                                  (substr(ssg.puck.f@cell_metadata$maxCelltype.l2, 4, 100) %in% l1.l2$name.l2)]
plt5 <- mySpatPlot1(ssg.puck.f, select_cells = sel.cells, 
                    cell_color = 'maxCelltype.l2', show_legend = T, 
                    cell_color_code = l2.colors, other_cell_color = 'gray95') +
  theme_noAxis(T) 
plt5 <- plt5 + ggtitle(paste0(puck, ' - ', 'maxCelltype.l2')) + 
  theme(plot.title = element_text(face = 'bold', vjust = -3, size = 13))
print(plt5)

#Plot Interstitium
to.use <- c("l2.VSMC","l2.VSMC/P","l2.FIB","l2.M-FIB","l2.MYOF","l2.aFIB",
            "l2.EC-PTC","l2.EC-AEA","l2.EC-DVR","l2.EC-AVR","l2.EC-LYM",
            "l2.MAC-M2","l2.T","l2.NKT","l2.N","l2.B","l2.MDC","l2.ncMON","l2.cDC","l2.pDC","l2.MAST","l2.PL")
ssg.puck.f <- subsetGiotto(ssg.puck, ssg.puck@cell_ID[ssg.puck@cell_metadata$maxCelltype.l2 %in% to.use])

sel.cells <- ssg.puck.f@cell_ID[(ssg.puck.f@cell_metadata$has.l2.rctd) & 
                                  (ssg.puck.f@cell_metadata$maxWeight.l2 >= 30) & 
                                  (substr(ssg.puck.f@cell_metadata$maxCelltype.l2, 4, 100) %in% l1.l2$name.l2)]
plt5 <- mySpatPlot1(ssg.puck.f, select_cells = sel.cells, 
                    cell_color = 'maxCelltype.l2', show_legend = T, 
                    cell_color_code = l2.colors, other_cell_color = 'gray95') +
  theme_noAxis(T) 
plt5 <- plt5 + ggtitle(paste0(puck, ' - ', 'maxCelltype.l2')) + 
  theme(plot.title = element_text(face = 'bold', vjust = -3, size = 13))
print(plt5)






###Fibrotic region 1: Plot Fig. 5d and ED Fig. 11a
ssg.puck <- ssg.list[[1]]
ssg.puck <- normalizeGiotto(ssg.puck, scalefactor = 10000, log_norm = T, scale_genes = T, scale_cells = F)
puck = "200903_01"

#Remove outlier beads
spatPlot(ssg.puck)
ssg.puck <- subsetGiottoLocs(ssg.puck, x_min=800, x_max=5700, y_min=300, y_max=5100)

#Plot fibrotic area all subclass.l2 cell types (ED Fig 13a)
spatPlot(ssg.puck)
ssg.puck.f <- subsetGiottoLocs(ssg.puck, x_min=500, x_max=2800, y_min=2300, y_max=4000)

sel.cells <- ssg.puck.f@cell_ID[(ssg.puck.f@cell_metadata$has.l2.rctd) & 
                                  (ssg.puck.f@cell_metadata$maxWeight.l2 >= 30) & 
                                  (substr(ssg.puck.f@cell_metadata$maxCelltype.l2, 4, 100) %in% l1.l2$name.l2)]
plt5 <- mySpatPlot1(ssg.puck.f, select_cells = sel.cells, 
                    cell_color = 'maxCelltype.l2', show_legend = T,
                    point_size = 2,
                    cell_color_code = l2.colors, other_cell_color = 'gray95') +
  theme_noAxis(T) 
plt5 <- plt5 + ggtitle(paste0(puck, ' - ', 'maxCelltype.l2')) + 
  theme(plot.title = element_text(face = 'bold', vjust = -3, size = 13))
print(plt5)

#Plot fibrotic area subclass.l2 main cell types (Fig. 5d)
unique(ssg.puck@cell_metadata$maxCelltype.l2)
to.use <- c("l2.POD","l2.EC-GC","l2.VSMC","l2.FIB","l2.EC-AEA","l2.MYOF","l2.cycMYOF",
            "l2.aFIB", "l2.PL","l2.MAC-M2","l2.MDC","l2.ncMON","l2.cycMNP","l2.NKT","l2.T","l2.cDC",
            "l2.aTAL1","l2.C-TAL")
ssg.puck.f <- subsetGiotto(ssg.puck, ssg.puck@cell_ID[ssg.puck@cell_metadata$maxCelltype.l2 %in% to.use])
ssg.puck.f <- subsetGiottoLocs(ssg.puck.f, x_min=500, x_max=2800, y_min=2300, y_max=4000)

sel.cells <- ssg.puck.f@cell_ID[(ssg.puck.f@cell_metadata$has.l2.rctd) & 
                                  (ssg.puck.f@cell_metadata$maxWeight.l2 >= 30) & 
                                  (substr(ssg.puck.f@cell_metadata$maxCelltype.l2, 4, 100) %in% l1.l2$name.l2)]
plt5 <- mySpatPlot1(ssg.puck.f, select_cells = sel.cells, 
                    cell_color = 'maxCelltype.l2', show_legend = T,
                    point_size = 2,
                    cell_color_code = l2.colors, other_cell_color = 'gray95') +
  theme_noAxis(T) 
plt5 <- plt5 + ggtitle(paste0(puck, ' - ', 'maxCelltype.l2')) + 
  theme(plot.title = element_text(face = 'bold', vjust = -3, size = 13))
print(plt5)

#Plot aTAL genes: ITGB6, EGF, SLC12A1  (Fig. 5d)
plt <- plotGenesRGB_DARK(ssg.puck.f, genes = c('ITGB6', 'EGF', 'SLC12A1'), min.expr = NULL,
                         pt.size = 2,
                         slot='scaled',
                         legend="right")
plt <- plt + ggtitle('ITGB6-EGF-SLC12A1') + theme(plot.title = element_text(face = 'bold', vjust = -3, size = 13))
print(plt)

#Plot MYOF/IMM cell types (Fig. 5d)
plt <- plotMetadata_DARK(ssg.puck.f, metaCols = c("l2.MYOF","l2.aTAL1","l2.T"), 
                         max.expr = NULL, min.expr = c(10,10,10), pt.size = 2,
                         hideAxis = F, legend = "right")
plt <- plt + ggtitle('MYOF-aTAL1-T') + theme(plot.title = element_text(face = 'bold', vjust = -3, size = 13))
print(plt)

#Plot marker genes: CD14, COL1A1, MYH11 (Fig. 5d)
plt <- plotGenesRGB_DARK(ssg.puck.f, genes = c('CD14', 'COL1A1', 'MYH11'), min.expr = NULL,
                         pt.size = 2,
                         slot='scaled',
                         legend="right") 
plt <- plt + ggtitle('CD14-COL1A1-MYH11') + theme(plot.title = element_text(face = 'bold', vjust = -3, size = 13))
print(plt)






###Fibrotic region 2: ED Fig. 11a
ssg.puck <- ssg.list[[3]]
ssg.puck <- normalizeGiotto(ssg.puck, scalefactor = 10000, log_norm = T, scale_genes = T, scale_cells = F)
puck = "210119_09"

#Remove outlier beads
spatPlot(ssg.puck)
ssg.puck <- subsetGiottoLocs(ssg.puck, x_min=500, x_max=5000, y_min=0, y_max=5000)

#Plot fibrotic area for subclass.l2 cell types (ED Fig 13a)
unique(ssg.puck@cell_metadata$maxCelltype.l2)
to.use <- unique(ssg.puck@cell_metadata$maxCelltype.l2)[!unique(ssg.puck@cell_metadata$maxCelltype.l2) %in% "l2.M-TAL"]
ssg.puck.f <- subsetGiotto(ssg.puck, ssg.puck@cell_ID[ssg.puck@cell_metadata$maxCelltype.l2 %in% to.use])
ssg.puck.f <- subsetGiottoLocs(ssg.puck.f, x_min=1500, x_max=3500, y_min=2200, y_max=3500)

plt5 <- mySpatPlot1(ssg.puck.f, select_cells = sel.cells, 
                    cell_color = 'maxCelltype.l2', show_legend = T,
                    point_size = 2,
                    cell_color_code = l2.colors, other_cell_color = 'gray95') +
  theme_noAxis(T) 
plt5 <- plt5 + ggtitle(paste0(puck, ' - ', 'maxCelltype.l2')) + 
  theme(plot.title = element_text(face = 'bold', vjust = -3, size = 13))
print(plt5)






###Fibrotic region 3: ED Fig. 11a
ssg.puck <- ssg.list[[4]]
ssg.puck <- normalizeGiotto(ssg.puck, scalefactor = 10000, log_norm = T, scale_genes = T, scale_cells = F)
puck = "210119_18"

#Remove outlier beads
spatPlot(ssg.puck)
ssg.puck <- subsetGiottoLocs(ssg.puck, x_min=100, x_max=4000, y_min=0, y_max=4100)

#Plot fibrotic area for subclass.l2 cell types (ED Fig 13a)
to.use <- unique(ssg.puck@cell_metadata$maxCelltype.l2)
ssg.puck.f <- subsetGiotto(ssg.puck, ssg.puck@cell_ID[ssg.puck@cell_metadata$maxCelltype.l2 %in% to.use])
ssg.puck.f <- subsetGiottoLocs(ssg.puck.f, x_min=1400, x_max=3500, y_min=3000, y_max=4100)
sel.cells <- ssg.puck.f@cell_ID[(ssg.puck.f@cell_metadata$has.l2.rctd) & 
                                  (ssg.puck.f@cell_metadata$maxWeight.l2 >= 30) & 
                                  (substr(ssg.puck.f@cell_metadata$maxCelltype.l2, 4, 100) %in% l1.l2$name.l2)]
plt5 <- mySpatPlot1(ssg.puck.f, select_cells = sel.cells, 
                    cell_color = 'maxCelltype.l2', show_legend = T,
                    point_size = 2,
                    cell_color_code = l2.colors, other_cell_color = 'gray95') +
  theme_noAxis(T) 
plt5 <- plt5 + ggtitle(paste0(puck, ' - ', 'maxCelltype.l2')) + 
  theme(plot.title = element_text(face = 'bold', vjust = -3, size = 13))
print(plt5)







###Combined dotplot of fibrotic areas (Fig 5d)
#"Puck_200903_01", "Puck_210119_09", "Puck_210119_18"
ssg.puck <- ssg.list[[1]]
ssg.puck.f1 <- subsetGiottoLocs(ssg.puck, x_min=500, x_max=2800, y_min=2300, y_max=4000)

ssg.puck <- ssg.list[[3]]
ssg.puck.f2 <- subsetGiottoLocs(ssg.puck, x_min=1500, x_max=3500, y_min=2200, y_max=3500)

ssg.puck <- ssg.list[[4]]
ssg.puck.f3 <- subsetGiottoLocs(ssg.puck, x_min=1400, x_max=3500, y_min=3000, y_max=4100)

ssg.list.f <- list(ssg.puck.f1,ssg.puck.f2,ssg.puck.f3)


counts.all <- t(merge.sparse(lapply(lapply(ssg.list.f, function(ss) ss@raw_exprs), t)))
meta.all <- do.call(rbind, lapply(ssg.list.f, function(ss) ss@cell_metadata))
ssg.all <- createGiottoObject(raw_exprs = counts.all, cell_metadata = meta.all, instructions = ginstuc)
ssg.all <- subsetGiotto(ssg.all, ssg.all@cell_ID[ssg.all@cell_metadata$inBoundary])

#Plot aTAL markers
genes2plot <- c("ESRRG","EGF","TMEM52B","KCTD1","TFAP2B","ITGB6","BICC1","CST3","TSPAN1","NRP1","TPM1","ERBB2",
                "MET","ARHGAP29","MACC1","PROM1")

rownames(ssg.all@cell_metadata) <- ssg.all@cell_metadata$cell_ID
sss <- CreateSeuratObject(counts = ssg.all@raw_exprs, meta.data = ssg.all@cell_metadata)
sss$maxCelltype.l2 <- substr(sss$maxCelltype.l2, 4, 100)

to.use <- c("aTAL1","C-TAL")
sss <- subset(sss, cells = rownames(sss@meta.data[sss@meta.data$maxCelltype.l2 %in% to.use,]))
sss$maxCelltype.l2 <- factor(sss$maxCelltype.l2, levels = c("aTAL1","C-TAL"))
sss <- SetIdent(sss, value = sss$maxCelltype.l2)

DotPlot(sss, features = (genes2plot), cols = c('gray95', 'orange'), ) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5, size = 9.4))

#Plot aSTR markers
genes2plot <- c("NEGR1","FBLN5","DCN","PDGFRA","LAMA2","FLRT2","LIFR","FGF14","COL16A1","COL14A1",
                "ACTA2","COL1A2","SORBS1",
                "COL1A1")

sss <- CreateSeuratObject(counts = ssg.all@raw_exprs, meta.data = ssg.all@cell_metadata)
sss$maxCelltype.l2 <- substr(sss$maxCelltype.l2, 4, 100)

to.use <- c("FIB","aFIB","MYOF")
sss <- subset(sss, cells = rownames(sss@meta.data[sss@meta.data$maxCelltype.l2 %in% to.use,]))
sss$maxCelltype.l2 <- factor(sss$maxCelltype.l2, levels = c("FIB","aFIB","MYOF"))
sss <- SetIdent(sss, value = sss$maxCelltype.l2)

DotPlot(sss, features = (genes2plot), cols = c('gray95', 'orange'), ) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5, size = 9.4))

#Plot IMM markers
genes2plot <- c("IGKC", "MZB1",
                "IL7R","CD247","CD96", 
                "STAB1", "F13A1", "MRC1", "CD163", 
                "HLA-DRA","CD14","SLC1A3")

sss <- CreateSeuratObject(counts = ssg.all@raw_exprs, meta.data = ssg.all@cell_metadata)
sss$maxCelltype.l2 <- substr(sss$maxCelltype.l2, 4, 100)

to.use <- c("PL","T","MAC-M2","MDC")
sss <- subset(sss, cells = rownames(sss@meta.data[sss@meta.data$maxCelltype.l2 %in% to.use,]))
sss$maxCelltype.l2 <- factor(sss$maxCelltype.l2, levels = c("PL","T","NKT","MAC-M2","MDC","ncMON"))
sss <- SetIdent(sss, value = sss$maxCelltype.l2)

DotPlot(sss, features = (genes2plot), cols = c('gray95', 'orange'), ) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5, size = 9.4))

#Barplot of Imm cell types
order <- c("l2.PL","l2.T","l2.NKT","l2.MAST","l2.MAC-M2","l2.MDC",
           "l2.ncMON","l2.cDC","l2.pDC")
prop1 <- (table(meta.all$maxCelltype.l2))
barplot(prop1[order],main = "Niche Cell Type Counts", cex.names = 0.5, 
        names.arg = order,las=2, col = as.character(l2.colors[order]))




###Full Tissue Puck - 200903_06 (ED Fig 4c)
ssg.puck <- ssg.list[[2]]
ssg.puck <- normalizeGiotto(ssg.puck, scalefactor = 10000, log_norm = T, scale_genes = T, scale_cells = F)
puck = "200903_06"

#Remove outlier beads
spatPlot(ssg.puck)
ssg.puck <- subsetGiottoLocs(ssg.puck, x_min=800, x_max=6000, y_min=0, y_max=5900)


#Plot subclass.l2 (ED Fig 4c)
sel.cells <- ssg.puck@cell_ID[(ssg.puck@cell_metadata$has.l2.rctd) & 
                                (ssg.puck@cell_metadata$maxWeight.l2 >= 30) & 
                                (substr(ssg.puck@cell_metadata$maxCelltype.l2, 4, 100) %in% l1.l2$name.l2)]
plt3 <- mySpatPlot1(ssg.puck, select_cells = sel.cells, 
                    cell_color = 'maxCelltype.l2', show_legend = T, 
                    cell_color_code = l2.colors, other_cell_color = 'gray95') +
  theme_noAxis(T) 
plt3 <- plt3 + ggtitle(paste0(puck, ' - ', 'maxCelltype.l2')) + 
  theme(plot.title = element_text(face = 'bold', vjust = -3, size = 13))
print(plt3)

#Plot Glom (ED Fig 4c)
to.use <- c("l2.POD","l2.PEC","l2.MC","l2.REN","l2.EC-AEA","l2.EC-GC")
ssg.puck.f <- subsetGiotto(ssg.puck, ssg.puck@cell_ID[ssg.puck@cell_metadata$maxCelltype.l2 %in% to.use])

sel.cells <- ssg.puck.f@cell_ID[(ssg.puck.f@cell_metadata$has.l2.rctd) & 
                                  (ssg.puck.f@cell_metadata$maxWeight.l2 >= 30) & 
                                  (substr(ssg.puck.f@cell_metadata$maxCelltype.l2, 4, 100) %in% l1.l2$name.l2)]
plt5 <- mySpatPlot1(ssg.puck.f, select_cells = sel.cells, 
                    cell_color = 'maxCelltype.l2', show_legend = T, 
                    cell_color_code = l2.colors, other_cell_color = 'gray95') +
  theme_noAxis(T) 
plt5 <- plt5 + ggtitle(paste0(puck, ' - ', 'maxCelltype.l2')) + 
  theme(plot.title = element_text(face = 'bold', vjust = -3, size = 13))
print(plt5)

#Plot tubules (ED Fig 4c)
to.use <- c("l2.CNT","l2.M-TAL","l2.M-PC","l2.IMCD","l2.DTL3","l2.DTL1","l2.ATL","l2.DTL2",
            "l2.C-PC","l2.aPT","l2.PT-S3","l2.DCT","l2.aTAL1","l2.C-TAL","l2.MD",
            "l2.PT-S1", "l2.IC-B", "l2.PT-S2","l2.PT-S3","l2.C-IC-A","l2.M-IC-A")
ssg.puck.f <- subsetGiotto(ssg.puck, ssg.puck@cell_ID[ssg.puck@cell_metadata$maxCelltype.l2 %in% to.use])

sel.cells <- ssg.puck.f@cell_ID[(ssg.puck.f@cell_metadata$has.l2.rctd) & 
                                  (ssg.puck.f@cell_metadata$maxWeight.l2 >= 30) & 
                                  (substr(ssg.puck.f@cell_metadata$maxCelltype.l2, 4, 100) %in% l1.l2$name.l2)]
plt5 <- mySpatPlot1(ssg.puck.f, select_cells = sel.cells, 
                    cell_color = 'maxCelltype.l2', show_legend = T, 
                    cell_color_code = l2.colors, other_cell_color = 'gray95') +
  theme_noAxis(T) 
plt5 <- plt5 + ggtitle(paste0(puck, ' - ', 'maxCelltype.l2')) + 
  theme(plot.title = element_text(face = 'bold', vjust = -3, size = 13))
print(plt5)

#Plot Interstitium (ED Fig 4c)
to.use <- c("l2.VSMC","l2.VSMC/P","l2.FIB","l2.M-FIB","l2.MYOF","l2.aFIB",
            "l2.EC-PTC","l2.EC-AEA","l2.EC-DVR","l2.EC-AVR","l2.EC-LYM",
            "l2.MAC-M2","l2.T","l2.NKT","l2.N","l2.B","l2.MDC","l2.ncMON","l2.cDC","l2.pDC","l2.MAST","l2.PL")
ssg.puck.f <- subsetGiotto(ssg.puck, ssg.puck@cell_ID[ssg.puck@cell_metadata$maxCelltype.l2 %in% to.use])

sel.cells <- ssg.puck.f@cell_ID[(ssg.puck.f@cell_metadata$has.l2.rctd) & 
                                  (ssg.puck.f@cell_metadata$maxWeight.l2 >= 30) & 
                                  (substr(ssg.puck.f@cell_metadata$maxCelltype.l2, 4, 100) %in% l1.l2$name.l2)]
plt5 <- mySpatPlot1(ssg.puck.f, select_cells = sel.cells, 
                    cell_color = 'maxCelltype.l2', show_legend = T, 
                    cell_color_code = l2.colors, other_cell_color = 'gray95') +
  theme_noAxis(T) 
plt5 <- plt5 + ggtitle(paste0(puck, ' - ', 'maxCelltype.l2')) + 
  theme(plot.title = element_text(face = 'bold', vjust = -3, size = 13))
print(plt5)



#Plot Renal Corpuscle region (Fig. 2e)
to.use <- c("l2.POD","l2.PEC","l2.MC","l2.REN","l2.MD","l2.EC-AEA","l2.EC-GC")
ssg.puck.f <- subsetGiotto(ssg.puck, ssg.puck@cell_ID[ssg.puck@cell_metadata$maxCelltype.l2 %in% to.use])
ssg.puck.f <- subsetGiottoLocs(ssg.puck.f, x_min=2100, x_max=3200, y_min=3000, y_max=3800)

plt5 <- mySpatPlot1(ssg.puck.f, 
                    cell_color = 'maxCelltype.l2', show_legend = T,
                    point_size = 2,
                    cell_color_code = l2.colors, other_cell_color = 'gray95') +
  theme_noAxis(T) 
plt5 <- plt5 + ggtitle(paste0(puck, ' - ', 'maxCelltype.l2')) + 
  theme(plot.title = element_text(face = 'bold', vjust = -3, size = 13))
print(plt5)

plt <- plotGenesRGB_DARK(ssg.puck.f, genes = c('REN', 'NOS1',"EMCN"), min.expr = NULL,
                         pt.size = 2,
                         slot='scaled',
                         legend="right")
plt <- plt + ggtitle('REN-NOS1-EMCN') + theme(plot.title = element_text(face = 'bold', vjust = -3, size = 13))
print(plt)

#Plot AEA subregion (ED Fig. 4f)
to.use <- c("l2.POD","l2.EC-GC","l2.VSMC","l2.FIB","l2.EC-AEA","l2.REN",
            "l2.C-TAL","l2.MD","l2.DCT","l2.C-PC",
            "l2.PT-S1","l2.PT-S2","l2.PT-S3")
ssg.puck.f <- subsetGiotto(ssg.puck, ssg.puck@cell_ID[ssg.puck@cell_metadata$maxCelltype.l2 %in% to.use])
ssg.puck.f <- subsetGiottoLocs(ssg.puck.f, x_min=4000, x_max=6000, y_min=2800, y_max=4000)

plt5 <- mySpatPlot1(ssg.puck.f, 
                    cell_color = 'maxCelltype.l2', show_legend = T,
                    point_size = 2,
                    cell_color_code = l2.colors, other_cell_color = 'gray95') +
  theme_noAxis(T) 
plt5 <- plt5 + ggtitle(paste0(puck, ' - ', 'maxCelltype.l2')) + 
  theme(plot.title = element_text(face = 'bold', vjust = -3, size = 13))
print(plt5)

plt <- plotGenesRGB_DARK(ssg.puck.f, genes = c('MYH11', 'SLC5A12',"NPHS2"), min.expr = NULL,
                         pt.size = 2,
                         slot='scaled',
                         legend="right")
plt <- plt + ggtitle('MYH11-SLC5A12-NPHS2') + theme(plot.title = element_text(face = 'bold', vjust = -3, size = 13))
print(plt)






###Plot Medullary pucks (ED Fig 4c)
pucks.use <- c("Puck_210113_37","Puck_210113_28")
ssg.list <- foreach(puck=pucks.use, .final=function(x) setNames(x, pucks.use)) %dopar% {
  i <- if(puck %in% pucks[[1]]) 1 else 2
  ssg <- readRDS(paste0(gioObjFolders[[i]], puck, ".rds"))
  ssg <- subsetGiotto(ssg, cell_ids = ssg@cell_ID[ssg@cell_metadata$has.l1.rctd])
  ssg <- subsetGiotto(ssg, cell_ids = ssg@cell_ID[ssg@cell_metadata$inBoundary])
  ssg 
}

##Outer Medulla - Puck_210113_37
ssg.puck <- ssg.list[[1]]
ssg.puck <- normalizeGiotto(ssg.puck, scalefactor = 10000, log_norm = T, scale_genes = T, scale_cells = F)
puck = "210113_37"

#Remove outlier beads
spatPlot(ssg.puck)
ssg.puck <- subsetGiottoLocs(ssg.puck, x_min=400, x_max=5000, y_min=200, y_max=5000)

#Plot subclass.l2 (ED Fig 4c)
sel.cells <- ssg.puck@cell_ID[(ssg.puck@cell_metadata$has.l2.rctd) & 
                                (ssg.puck@cell_metadata$maxWeight.l2 >= 30) & 
                                (substr(ssg.puck@cell_metadata$maxCelltype.l2, 4, 100) %in% l1.l2$name.l2)]
plt3 <- mySpatPlot1(ssg.puck, select_cells = sel.cells, 
                    cell_color = 'maxCelltype.l2', show_legend = T, 
                    cell_color_code = l2.colors, other_cell_color = 'gray95') +
  theme_noAxis(T) 
plt3 <- plt3 + ggtitle(paste0(puck, ' - ', 'maxCelltype.l2')) + 
  theme(plot.title = element_text(face = 'bold', vjust = -3, size = 13))
print(plt3)

#Plot Gloms (ED Fig 4c)
to.use <- c("l2.POD","l2.PEC","l2.MC","l2.REN","l2.EC-AEA","l2.EC-GC")
ssg.puck.f <- subsetGiotto(ssg.puck, ssg.puck@cell_ID[ssg.puck@cell_metadata$maxCelltype.l2 %in% to.use])

sel.cells <- ssg.puck.f@cell_ID[(ssg.puck.f@cell_metadata$has.l2.rctd) & 
                                  (ssg.puck.f@cell_metadata$maxWeight.l2 >= 30) & 
                                  (substr(ssg.puck.f@cell_metadata$maxCelltype.l2, 4, 100) %in% l1.l2$name.l2)]
plt5 <- mySpatPlot1(ssg.puck.f, select_cells = sel.cells, 
                    cell_color = 'maxCelltype.l2', show_legend = T, 
                    cell_color_code = l2.colors, other_cell_color = 'gray95') +
  theme_noAxis(T) 
plt5 <- plt5 + ggtitle(paste0(puck, ' - ', 'maxCelltype.l2')) + 
  theme(plot.title = element_text(face = 'bold', vjust = -3, size = 13))
print(plt5)

#Plot tubules (ED Fig 4c)
to.use <- c("l2.CNT","l2.M-TAL","l2.M-PC","l2.IMCD","l2.DTL3","l2.DTL1","l2.ATL","l2.DTL2",
            "l2.C-PC","l2.aPT","l2.PT-S3","l2.DCT","l2.aTAL1","l2.C-TAL","l2.MD",
            "l2.PT-S1", "l2.IC-B", "l2.PT-S2","l2.PT-S3","l2.C-IC-A","l2.M-IC-A")
ssg.puck.f <- subsetGiotto(ssg.puck, ssg.puck@cell_ID[ssg.puck@cell_metadata$maxCelltype.l2 %in% to.use])

sel.cells <- ssg.puck.f@cell_ID[(ssg.puck.f@cell_metadata$has.l2.rctd) & 
                                  (ssg.puck.f@cell_metadata$maxWeight.l2 >= 30) & 
                                  (substr(ssg.puck.f@cell_metadata$maxCelltype.l2, 4, 100) %in% l1.l2$name.l2)]
plt5 <- mySpatPlot1(ssg.puck.f, select_cells = sel.cells, 
                    cell_color = 'maxCelltype.l2', show_legend = T, 
                    cell_color_code = l2.colors, other_cell_color = 'gray95') +
  theme_noAxis(T) 
plt5 <- plt5 + ggtitle(paste0(puck, ' - ', 'maxCelltype.l2')) + 
  theme(plot.title = element_text(face = 'bold', vjust = -3, size = 13))
print(plt5)

#Plot Interstitium
to.use <- c("l2.VSMC","l2.VSMC/P","l2.FIB","l2.M-FIB","l2.MYOF","l2.aFIB",
            "l2.EC-PTC","l2.EC-AEA","l2.EC-DVR","l2.EC-AVR","l2.EC-LYM",
            "l2.MAC-M2","l2.T","l2.NKT","l2.N","l2.B","l2.MDC","l2.ncMON","l2.cDC","l2.pDC","l2.MAST","l2.PL")
ssg.puck.f <- subsetGiotto(ssg.puck, ssg.puck@cell_ID[ssg.puck@cell_metadata$maxCelltype.l2 %in% to.use])

sel.cells <- ssg.puck.f@cell_ID[(ssg.puck.f@cell_metadata$has.l2.rctd) & 
                                  (ssg.puck.f@cell_metadata$maxWeight.l2 >= 30) & 
                                  (substr(ssg.puck.f@cell_metadata$maxCelltype.l2, 4, 100) %in% l1.l2$name.l2)]
plt5 <- mySpatPlot1(ssg.puck.f, select_cells = sel.cells, 
                    cell_color = 'maxCelltype.l2', show_legend = T, 
                    cell_color_code = l2.colors, other_cell_color = 'gray95') +
  theme_noAxis(T) 
plt5 <- plt5 + ggtitle(paste0(puck, ' - ', 'maxCelltype.l2')) + 
  theme(plot.title = element_text(face = 'bold', vjust = -3, size = 13))
print(plt5)



##Inner Medulla - Puck_210113_28
ssg.puck <- ssg.list[[2]]
ssg.puck <- normalizeGiotto(ssg.puck, scalefactor = 10000, log_norm = T, scale_genes = T, scale_cells = F)
puck = "210113_28"

#Remove outlier beads
spatPlot(ssg.puck)
ssg.puck <- subsetGiottoLocs(ssg.puck, x_min=300, x_max=5000, y_min=0, y_max=4100)

#Plot subclass.l2 (ED Fig 4c)
sel.cells <- ssg.puck@cell_ID[(ssg.puck@cell_metadata$has.l2.rctd) & 
                                (ssg.puck@cell_metadata$maxWeight.l2 >= 30) & 
                                (substr(ssg.puck@cell_metadata$maxCelltype.l2, 4, 100) %in% l1.l2$name.l2)]
plt3 <- mySpatPlot1(ssg.puck, select_cells = sel.cells, 
                    cell_color = 'maxCelltype.l2', show_legend = T, 
                    cell_color_code = l2.colors, other_cell_color = 'gray95') +
  theme_noAxis(T) 
plt3 <- plt3 + ggtitle(paste0(puck, ' - ', 'maxCelltype.l2')) + 
  theme(plot.title = element_text(face = 'bold', vjust = -3, size = 13))
print(plt3)

#Plot Gloms (ED Fig 4c)
to.use <- c("l2.POD","l2.PEC","l2.MC","l2.REN","l2.EC-AEA","l2.EC-GC")
ssg.puck.f <- subsetGiotto(ssg.puck, ssg.puck@cell_ID[ssg.puck@cell_metadata$maxCelltype.l2 %in% to.use])

sel.cells <- ssg.puck.f@cell_ID[(ssg.puck.f@cell_metadata$has.l2.rctd) & 
                                  (ssg.puck.f@cell_metadata$maxWeight.l2 >= 30) & 
                                  (substr(ssg.puck.f@cell_metadata$maxCelltype.l2, 4, 100) %in% l1.l2$name.l2)]
plt5 <- mySpatPlot1(ssg.puck.f, select_cells = sel.cells, 
                    cell_color = 'maxCelltype.l2', show_legend = T, 
                    cell_color_code = l2.colors, other_cell_color = 'gray95') +
  theme_noAxis(T) 
plt5 <- plt5 + ggtitle(paste0(puck, ' - ', 'maxCelltype.l2')) + 
  theme(plot.title = element_text(face = 'bold', vjust = -3, size = 13))
print(plt5)

#Plot tubules (ED Fig 4c)
to.use <- c("l2.CNT","l2.M-TAL","l2.M-PC","l2.IMCD","l2.DTL3","l2.DTL1","l2.ATL","l2.DTL2",
            "l2.C-PC","l2.aPT","l2.PT-S3","l2.DCT","l2.aTAL1","l2.C-TAL","l2.MD",
            "l2.PT-S1", "l2.IC-B", "l2.PT-S2","l2.PT-S3","l2.C-IC-A","l2.M-IC-A")
ssg.puck.f <- subsetGiotto(ssg.puck, ssg.puck@cell_ID[ssg.puck@cell_metadata$maxCelltype.l2 %in% to.use])

sel.cells <- ssg.puck.f@cell_ID[(ssg.puck.f@cell_metadata$has.l2.rctd) & 
                                  (ssg.puck.f@cell_metadata$maxWeight.l2 >= 30) & 
                                  (substr(ssg.puck.f@cell_metadata$maxCelltype.l2, 4, 100) %in% l1.l2$name.l2)]
plt5 <- mySpatPlot1(ssg.puck.f, select_cells = sel.cells, 
                    cell_color = 'maxCelltype.l2', show_legend = T, 
                    cell_color_code = l2.colors, other_cell_color = 'gray95') +
  theme_noAxis(T) 
plt5 <- plt5 + ggtitle(paste0(puck, ' - ', 'maxCelltype.l2')) + 
  theme(plot.title = element_text(face = 'bold', vjust = -3, size = 13))
print(plt5)

#Plot Interstitium (ED Fig 4c)
to.use <- c("l2.VSMC","l2.VSMC/P","l2.FIB","l2.M-FIB","l2.MYOF","l2.aFIB",
            "l2.EC-PTC","l2.EC-AEA","l2.EC-DVR","l2.EC-AVR","l2.EC-LYM",
            "l2.MAC-M2","l2.T","l2.NKT","l2.N","l2.B","l2.MDC","l2.ncMON","l2.cDC","l2.pDC","l2.MAST","l2.PL")
ssg.puck.f <- subsetGiotto(ssg.puck, ssg.puck@cell_ID[ssg.puck@cell_metadata$maxCelltype.l2 %in% to.use])

sel.cells <- ssg.puck.f@cell_ID[(ssg.puck.f@cell_metadata$has.l2.rctd) & 
                                  (ssg.puck.f@cell_metadata$maxWeight.l2 >= 30) & 
                                  (substr(ssg.puck.f@cell_metadata$maxCelltype.l2, 4, 100) %in% l1.l2$name.l2)]
plt5 <- mySpatPlot1(ssg.puck.f, select_cells = sel.cells, 
                    cell_color = 'maxCelltype.l2', show_legend = T, 
                    cell_color_code = l2.colors, other_cell_color = 'gray95') +
  theme_noAxis(T) 
plt5 <- plt5 + ggtitle(paste0(puck, ' - ', 'maxCelltype.l2')) + 
  theme(plot.title = element_text(face = 'bold', vjust = -3, size = 13))
print(plt5)


###
