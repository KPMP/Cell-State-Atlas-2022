library(doParallel)
library(foreach)
library(Seurat)
library(dplyr)
library(tidyr)
registerDoParallel(5)
source("./utils.R")


plot.dir <- "./Plots/"

gioObjFolder <- paste0("./deconv_results/GiottoObjects/")
l1.l2 <- read.table(paste0("./deconv_results/HuBMAP_l1-l2_subclasses.tsv"), header = T)

# prepping puck names and addressed
pucks <- list.files(gioObjFolder)
pucks <- pucks[grepl("Puck_[0-9]+_\\d+\\.rds", pucks)]
pucks <- substr(pucks, 1, nchar(pucks) - 4)

# Giotto's stuff
python_path <- "/media/Home_Raid1_Voyager/kian/anaconda3/envs/giotto_200703/bin/python"
ginstuc <- createGiottoInstructions(python_path = python_path,
                                    show_plot = F, save_plot = F)
# Reading the region information
region_info <- readxl::read_xlsx("./SlideSeq_puck_regions.xlsx", )

# setting up the color table
colorTable <- data.table::fread("./Cluster_Color_Table.tsv", select = 1:20)
l2.colors.df <- colorTable[, c('subclass.l2_label', 'subclass.l2_color')]
l2.colors.df <- l2.colors.df[!duplicated(l2.colors.df), ]
l2.colors.df <- l2.colors.df[l2.colors.df$subclass.l2_label %in% l1.l2$name.l2]
l2.colors <- l2.colors.df$subclass.l2_color
names(l2.colors) <- paste0('l2.', l2.colors.df$subclass.l2_label)


# plotting panels in Extended Data Figure 5a: ####
plot_umi_weights_counts <- function(gObj.list, file, color.vec, titleReg){
  maxCts <- unlist(sapply(gObj.list, function(ssg) setNames(ssg@cell_metadata$maxCelltype.l2, ssg@cell_ID)))
  maxWs <- unlist(sapply(gObj.list, function(ssg) setNames(ssg@cell_metadata$maxWeight.l2, ssg@cell_ID)))
  nUMIs <- unlist(sapply(gObj.list, function(ssg) setNames(ssg@cell_metadata$nUMI, ssg@cell_ID)))
  
  df <- data.frame(maxCelltype = maxCts, maxWeight = maxWs, nUMI = nUMIs, row.names = names(nUMIs))
  df <- df[substr(df$maxCelltype, 4, 100) %in% l1.l2$name.l2, ]
  df$maxCelltype <- factor(df$maxCelltype, levels = names(color.vec))
  
  ctCounts <- table(df$maxCelltype)
  ct2choose <- names(ctCounts)[ctCounts >= 10]
  
  count.df <- data.frame(ctCounts[ct2choose])
  colnames(count.df) <- c('celltype', 'count')
  count.df$celltype <- factor(count.df$celltype, levels = names(color.vec))
  
  df2p <- df[df$maxCelltype %in% ct2choose, ]
  
  plt1 <- ggplot(df2p, aes(x = maxCelltype, y = nUMI, fill = maxCelltype)) + 
            geom_violin(draw_quantiles = c(0.5), scale = "width") + 
            ggtitle(paste0("UMI counts in " , titleReg)) + 
            theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5), legend.position = "None", 
                  plot.title = element_text(hjust=0.5, face='bold')) + 
            scale_fill_manual(values = color.vec) +
            scale_y_continuous(trans='log10', limits = c(100, 5000))
  
  
  plt2 <- ggplot(df2p, aes(x = maxCelltype, y = maxWeight, fill = maxCelltype)) + 
            geom_violin(draw_quantiles = c(0.5), scale = "width", ) + 
            ggtitle(paste0("Maximum weights in " , titleReg)) + 
            theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5), legend.position = "None", 
                  plot.title = element_text(hjust=0.5, face='bold')) + 
            scale_fill_manual(values = color.vec)
  
  plt3 <- ggplot(count.df, aes(x=celltype, y = count, fill=celltype)) +
            geom_bar(stat='identity') + 
            ggtitle(paste0("Cell type counts in " , titleReg)) + 
            theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5), legend.position = "None", 
                  plot.title = element_text(hjust=0.5, face='bold')) + 
            scale_fill_manual(values = color.vec) #+ 
            # scale_y_continuous(trans='log10')
  pdf(file, width = 10)
  print(plt1)
  print(plt2)
  print(plt3)
  dev.off()

}

puck.sel <- pucks
ssg.sel <- foreach(puck=puck.sel) %dopar%{
  ssg <- readRDS(paste0(gioObjFolder, puck, ".rds"))
  ssg <- filterWeights(ssg, 'maxWeight.l2', 'maxCelltype.l2', 30, c() )
  ssg
}
# 1st, 2nd and 4th panels
plot_umi_weights_counts(ssg.sel, paste0(plot.dir, "ED5_a_1.pdf"), l2.colors,
                        titleReg='All Pucks')
# 3rd panel (normalized cell type distribution between regions)
getCtCounts <- function(gObj.list){
  maxCts <- unlist(sapply(gObj.list, function(ssg) setNames(ssg@cell_metadata$maxCelltype.l2, ssg@cell_ID)))
  maxCts <- maxCts[substr(maxCts, 4, 100) %in% l1.l2$name.l2]
  ctCounts <- table(maxCts)
  return(ctCounts)  
}

reg <- 'Cortex'
puck.sel <- intersect(region_info$Specimen[region_info$`Region Level 1` == reg], unlist(pucks))
ssg.sel <- foreach(puck=puck.sel) %dopar%{
  ssg <- readRDS(paste0(gioObjFolder, puck, ".rds"))
  ssg <- filterWeights(ssg, 'maxWeight.l2', 'maxCelltype.l2', 30, c() )
  ssg
}

cort.df <- getCtCounts(ssg.sel)
cort.df <- as.data.frame(cort.df)
cort.df$Freq <- cort.df$Freq / sum(cort.df$Freq) * 100


reg <- 'Outer Medulla'
puck.sel <- intersect(region_info$Specimen[region_info$`Region Level 2` == reg], unlist(pucks))
ssg.sel <- foreach(puck=puck.sel) %dopar%{
  ssg <- readRDS(paste0(gioObjFolder, puck, ".rds"))
  ssg <- filterWeights(ssg, 'maxWeight.l2', 'maxCelltype.l2', 30, c() )
  ssg
}

OM.df <- getCtCounts(ssg.sel)
OM.df <- as.data.frame(OM.df)
OM.df$Freq <- OM.df$Freq / sum(OM.df$Freq) * 100


reg <- 'Inner Medulla'
puck.sel <- intersect(region_info$Specimen[region_info$`Region Level 2` == reg], unlist(pucks))
ssg.sel <- foreach(puck=puck.sel) %dopar%{
  ssg <- readRDS(paste0(gioObjFolder, puck, ".rds"))
  ssg <- filterWeights(ssg, 'maxWeight.l2', 'maxCelltype.l2', 30, c() )
  ssg
}

IM.df <- getCtCounts(ssg.sel)
IM.df <- as.data.frame(IM.df)
IM.df$Freq <- IM.df$Freq / sum(IM.df$Freq) * 100

prop.df <- dplyr::full_join(cort.df, OM.df, by="maxCts")

prop.df <- Reduce(function(x, y) dplyr::full_join(x, y, by="maxCts"), list(cort.df, OM.df, IM.df))
colnames(prop.df)[2:4] <- c("Cortex", "Outer Medulla", "Inner Medulla")
prop.df$maxCts <- factor(prop.df$maxCts, levels = names(l2.colors))

prop.df.norm <- prop.df
prop.df.norm[, c('Cortex', 'Outer Medulla', 'Inner Medulla')] <- 100 * sweep(prop.df.norm[, c('Cortex', 'Outer Medulla', 'Inner Medulla')], 1, rowSums(prop.df.norm[, c('Cortex', 'Outer Medulla', 'Inner Medulla')]), "/")

pdf(paste0(plot.dir, "ED5_a_2.pdf"), width = 9, height=5)
ggplot(reshape2::melt(prop.df.norm), aes(x = maxCts, y = value, fill = variable)) + 
  geom_bar(stat='identity', width = 0.8) + 
  scale_fill_manual(values = c(Cortex="#508F5B", `Outer Medulla`="#3D8ACC", `Inner Medulla`="#7014CC")) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5))
dev.off()

# Extended Data Figure 5b : marker gene dot plot for slide-seq ####
genes2plot.df <- read.table("./Markers_dotplot.txt", sep = "\t")
l2classes2 <- genes2plot.df$V1

ssg.list <- foreach(puck=unlist(pucks))%dopar%{
  ssg <- readRDS(paste0(gioObjFolder, puck, ".rds"))
  ssg <- filterWeights(ssg, 'maxWeight.l2', 'maxCelltype.l2', 30, c() )
  ssg 
}

counts.all <- t(merge.sparse(lapply(lapply(ssg.list, function(ss) ss@raw_exprs), t)))
meta.all <- do.call(rbind, lapply(ssg.list, function(ss) ss@cell_metadata))
ssg.all <- createGiottoObject(raw_exprs = counts.all, cell_metadata = meta.all, instructions = ginstuc)
ssg.all <- filterWeights(ssg.all, 'maxWeight.l2', 'maxCelltype.l2', 30, c())

genes2plot <- unlist(lapply(genes2plot.df$V2, function(x) strsplit(x, ", ")))
genes2plot <- unique(genes2plot[genes2plot %in% ssg.all@gene_ID])

rownames(ssg.all@cell_metadata) <- ssg.all@cell_metadata$cell_ID
sss.norm <- CreateSeuratObject(counts = ssg.all@raw_exprs, meta.data = ssg.all@cell_metadata)
sss.norm <- NormalizeData(sss.norm, normalization.method = 'LogNormalize')
sss.norm$maxCelltype.l2 <- substr(sss.norm$maxCelltype.l2, 4, 100)

ct.toremove <- unique(sss.norm$maxCelltype.l2[!sss.norm$maxCelltype.l2 %in% l2classes2])
sss.norm <- subset(sss.norm, cells = colnames(sss.norm)[!sss.norm$maxCelltype.l2 %in% ct.toremove])
sss.norm$maxCelltype.l2 <- factor(sss.norm$maxCelltype.l2, levels = l2classes2)
sss.norm <- SetIdent(sss.norm, value = sss.norm$maxCelltype.l2)

pdf(paste0(plot.dir, "ED5_b.pdf"), height=8, width = 16)
DotPlot(sss.norm, features = (genes2plot), cols = c('gray95', 'orange'), scale = F) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5, size = 9.4))
dev.off()

# Figure 2h: Renal corpuscle dot plot####
genes2plot.glom <- c("NPHS2","EMCN","POSTN","VCAM1","SLC5A12","PALMD","MYH11","REN","NOS1")
ct2plot <- c("POD","EC-GC","MC","PEC","PT-S1/2","EC-AEA","VSMC","REN","MD")

sel.pucks <- region_info$Specimen[region_info$`Region Level 1` == 'Cortex']
ssg.list <- foreach(puck=sel.pucks)%dopar%{
  print(puck)
  ssg <- readRDS(paste0(gioObjFolder, puck, ".rds"))
  ssg <- filterWeights(ssg, 'maxWeight.l2', 'maxCelltype.l2', 30, c() )
  ssg 
}

counts.all <- t(merge.sparse(lapply(lapply(ssg.list, function(ss) ss@raw_exprs), t)))
meta.all <- do.call(rbind, lapply(ssg.list, function(ss) ss@cell_metadata))
ssg.all <- createGiottoObject(raw_exprs = counts.all, cell_metadata = meta.all, instructions = ginstuc)
ssg.all <- filterWeights(ssg.all, 'maxWeight.l2', 'maxCelltype.l2', 30, c())

rownames(ssg.all@cell_metadata) <- ssg.all@cell_metadata$cell_ID
sss <- CreateSeuratObject(counts = ssg.all@raw_exprs, meta.data = ssg.all@cell_metadata)
maxCt.ptmix <- sss$maxCelltype.l2
maxCt.ptmix <- substr(maxCt.ptmix, 4, 100)
maxCt.ptmix[maxCt.ptmix %in% c('PT-S1', 'PT-S2')] <- 'PT-S1/2'
sss[['maxCt.ptmix']] <- maxCt.ptmix
sss$maxCt.ptmix <- factor(sss$maxCt.ptmix, levels = ct2plot)
sss <- SetIdent(sss, value = 'maxCt.ptmix')
sss.sub <- NormalizeData(sss, normalization.method = 'LogNormalize')
sss.sub <- subset(sss.sub, idents = ct2plot)


pdf(paste0(plot.dir, "/Fig2_h.pdf"), height=5, width = 7)
DotPlot(sss.sub, features = (genes2plot.glom), cols = c('gray95', 'orange'), dot.scale = 15) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5, size = 9.4))
dev.off()

# Figure 2d: outer to inner medulla transition ####
rotateGiotto <- function(gObj, ang){
  R <- matrix(c(cos(ang), -sin(ang), sin(ang), cos(ang)), 2, byrow = T)
  gObj@spatial_locs[, c(1, 2)] <- as.data.frame(t(R %*% t(subset(gObj@spatial_locs, 
                                                                 select = c("sdimx", "sdimy")))))
  gObj
  
}
trans.puck <- c("Puck_210113_23")
ssg <- readRDS(paste0(gioObjFolder, trans.puck, ".rds"))
ssg <- subsetGiotto(ssg, cell_ids = ssg@cell_ID[ssg@cell_metadata$has.l1.rctd])
ssg <- subsetGiotto(ssg, cell_ids = ssg@cell_ID[ssg@cell_metadata$inBoundary])
ssg <- normalizeGiotto(ssg, scalefactor = 10000, log_norm = T, scale_genes = T, scale_cells = F)

ssg <- rotateGiotto(ssg, ang=-28 / 180 * pi)

plt3 <- plotMetadata_DARK(ssg, metaCols = c("l2.M-TAL", "l2.ATL", "l2.EC-DVR"),
                        min.expr = rep(20, 3), max.expr = rep(80, 3),
                        pt.size = 0.2) +
  plot_annotation(title=trans.puck,
                  theme = theme(plot.title = element_text(face='bold', hjust = 0.5)))
plt4 <- plotGenesRGB_DARK(ssg, genes = c('SLC12A1', 'SH3GL3', 'SLC14A1'), slot = 'scaled',
                      max.expr = c(2.5, 5, 4),
                      min.expr = c(0, 0, 0),
                      pt.size = 0.2) +
  plot_annotation(title=trans.puck,
                  theme = theme(plot.title = element_text(face='bold', hjust = 0.5)))
pdf(paste0(plot.dir, "Fig2d_IM-OM.pdf"), width = 10, height=10)
print(plt3)
print(plt4)
dev.off()


# Figure 2d: mean cell type proportion per region####
# reading all Giotto objects
maxW.l2.df <- foreach(puck = unlist(pucks), .final=function(x) setNames(x, unlist(pucks)))%dopar%{
  print(puck)
  ssg <- readRDS(paste0(gioObjFolder, puck, ".rds"))
  ssg <- filterWeights(ssg, 'maxWeight.l2', 'maxCelltype.l2', 30, c() )
  ssg <- subsetGiotto(ssg, cell_ids = ssg@cell_ID[ssg@cell_metadata$inBoundary])
  ssg <- subsetGiotto(ssg, cell_ids = ssg@cell_ID[substr(ssg@cell_metadata$maxCelltype.l2, 4, 100) %in% l1.l2$name.l2])
  ssg@cell_metadata[, c('cell_ID', 'maxCelltype.l2', 'maxWeight.l2')]
}

maxW.l2.df <- do.call(rbind, maxW.l2.df)
maxW.l2.df$puck.id <- sapply(strsplit(maxW.l2.df$cell_ID, '-'), function(x) x[1])
maxW.l2.df$maxCelltype.l2 <- factor(maxW.l2.df$maxCelltype.l2,
                            levels = names(l2.colors))

ct.puck.count <- maxW.l2.df %>% count(puck.id, maxCelltype.l2) # counting celltypes per puck
ct.puck.count <- ct.puck.count %>% group_by(puck.id) %>% mutate(inPuck.frac = n / sum(n)) # fraction of cell type per puck
ct.puck.count <- ct.puck.count %>% group_by(maxCelltype.l2) %>% mutate(norm.frac = inPuck.frac / max(inPuck.frac)) # normalizing fraction per puck across puck per cell type


ct.puck.count2 <- cbind(ct.puck.count, 
                       do.call(rbind, lapply(ct.puck.count$puck.id, function(pc){ 
                         region_info[region_info$Specimen==pc, c("Region Level 1", "Region Level 2", "Data Source")]
                         })))
region.order <- c("Cortex", "Outer Medulla", "Inner/Outer Medulla", "Inner Medulla")

regMean <- ct.puck.count2 %>% filter(`Region Level 2` != "Cortex / Medulla") %>% 
  group_by(maxCelltype.l2, `Region Level 2`) %>%
  summarise(regMeanFrac=mean(inPuck.frac))
regMean <- regMean %>% group_by(maxCelltype.l2) %>%
  mutate(normRegionMean = regMeanFrac / max(regMeanFrac))

regMean.wide <- pivot_wider(regMean[c('Region Level 2', 'maxCelltype.l2', 'normRegionMean')], 
                                   names_from = 'Region Level 2', values_from = 'normRegionMean')
regMean.wide[is.na(regMean.wide)] <- 0 # filling 0 for unobserved cell types

toPlot.df <- reshape2::melt(regMean.wide, variable.name = 'region', value.name = 'frac.norm')
toPlot.df$region <- factor(toPlot.df$region, levels = region.order)

# Selecting a subset of cell types 
sel.ctypes <- c('l2.PEC', 'l2.POD', 'l2.MC', 'l2.EC-GC', 'l2.REN', 'l2.EC-AEA', 
                'l2.VSMC', 'l2.PT-S1', 'l2.PT-S2', 'l2.PT-S3', 'l2.aPT', 
                'l2.EC-PTC', 'l2.EC-LYM', 'l2.FIB', 'l2.MYOF', 'l2.aFIB', 
                'l2.C-TAL', 'l2.aTAL1', 'l2.MD', 'l2.DCT', 'l2.CNT', 
                'l2.C-PC', 'l2.C-IC-A', 'l2.IC-B', 'l2.DTL1', 'l2.DTL2', 
                'l2.EC-DVR', 'l2.VSMC/P', 'l2.M-TAL', 'l2.M-PC', 'l2.M-IC-A', 
                'l2.EC-AVR', 'l2.DTL3', 'l2.ATL', 'l2.IMCD', 'l2.M-FIB' )

toPlot.df <- toPlot.df[toPlot.df$maxCelltype.l2 %in% sel.ctypes, ]
toPlot.df$maxCelltype.l2 <- factor(toPlot.df$maxCelltype.l2, levels = sel.ctypes)

pdf(paste0(plot.dir, "/", "Fig2d_celltypeProp.pdf"), width = 3.3, height = 8)
ggplot(data = toPlot.df, aes(x = region, y = maxCelltype.l2, fill = frac.norm)) + 
  geom_tile() + theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1)) + 
  scale_fill_continuous(type='viridis')
dev.off()
