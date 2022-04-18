# snCv3 - Degenerative State Gene Set and Score (Supplementary Table 12) -----------------------------------------
library(Seurat)
library(SeuratDisk)
library(dplyr)
library(Matrix)
library(ggplot2)
library(forcats)
library(cowplot)
library(patchwork)

KBR <- LoadH5Seurat("Kidney_Healthy-Injury_Cell_Atlas_snCv3_Seurat_03282022.h5Seurat")
Idents(object = KBR) <- "clusters"
Idents(object = KBR) <- factor(Idents(object = KBR), levels = 1:119)
KBR.sub <- subset(KBR, idents = 1:100)

load("color_factors.robj")


#Load Function
getconservedDEGs <- function(seurat.obj, celltype) {
  cons.DEGs <- do.call(rbind, lapply(celltype, function(ct) {
    print(paste("Running for cell type:", ct))
    
    mark <- FindConservedMarkers(seurat.obj, ident.1 = ct, grouping.var = "condition.l1",
                                 min.pct = 0.25, only.pos = TRUE, max.cells.per.ident = 300)
    if (length(rownames(mark))>0) {
      mark$gene <- rownames(mark)
      mark$cluster <- rep(ct, length(rownames(mark)))
      rownames(mark) <- paste(mark$cluster, mark$gene, sep = ".")
      return(mark)
    } else {
      print(paste("No DEGs for cell type:", ct))
    }
  }))
  return(cons.DEGs)
}


#Prepare Seurat object
Idents(KBR.sub) <- "subclass.l3"
Idents(KBR.sub) <- factor(Idents(KBR.sub), levels = c(
  "POD","dPOD","PEC","PT-S1","PT-S2","PT-S3","aPT","cycPT","dPT","dPT/DTL",
  "DTL2","DTL1","DTL3","dDTL3","ATL","dATL","aTAL1","aTAL2","M-TAL","dM-TAL",
  "C-TAL","dC-TAL","MD","DCT1","DCT2","dDCT","cycDCT","CNT","CNT-PC","dCNT",
  "cycCNT","CCD-PC","OMCD-PC","dOMCD-PC","IMCD","dIMCD","CCD-IC-A","CNT-IC-A",
  "dC-IC-A","OMCD-IC-A","tPC-IC","IC-B","PapE","EC-GC","EC-AEA","EC-DVR","EC-PTC",
  "dEC-PTC","EC-AVR","dEC","cycEC","EC-LYM","MC","REN","VSMC","VSMC/P","dVSMC",
  "MYOF","cycMYOF","FIB","M-FIB","dM-FIB","aFIB","dFIB","B","PL","T","NKT","MAST",
  "MAC-M2","cycMNP","MDC","cDC","pDC","ncMON","N","SC/NEU"
))
table(Idents(KBR.sub))





###Conserved Degen State Markers - POD
KBR.sub.pod <- subset(KBR.sub, idents = c("POD","dPOD"))

celltype <- levels(Idents(KBR.sub.pod))

cons.POD.DEGs <- getconservedDEGs(KBR.sub.pod, celltype)



###Conserved Degen State Markers - PT S1-S3
KBR.sub.pts13 <- subset(KBR.sub, idents = c("PT-S1","PT-S2","PT-S3",
                                            "dPT","dPT/DTL"))
KBR.sub.pts13 <- RenameIdents(object = KBR.sub.pts13, 
                              "PT-S1" = "PT", 
                              "PT-S2" = "PT",
                              "PT-S3" = "PT",
                              "dPT/DTL" = "dPT")
KBR.sub.pts13$subclass.l3 <- Idents(KBR.sub.pts13)

celltype <- levels(Idents(KBR.sub.pts13))
cons.pts13.DEGs <- getconservedDEGs(KBR.sub.pts13, celltype)



###Conserved Degen State Markers - DTL3 - ATL
KBR.sub.dtl3 <- subset(KBR.sub, idents = c("DTL3", "dDTL3", "ATL", "dATL"))
KBR.sub.dtl3 <- RenameIdents(object = KBR.sub.dtl3, 
                             "DTL3" = "DTL-ATL",
                             "dDTL3" = "dDTL-ATL",
                             "ATL" = "DTL-ATL",
                             "dATL" = "dDTL-ATL")

celltype <- levels(Idents(KBR.sub.dtl3))
cons.dtl3.DEGs <- getconservedDEGs(KBR.sub.dtl3, celltype)



###Conserved Degen State Markers - TAL
KBR.sub.tal <- subset(KBR.sub, idents = c("M-TAL","dM-TAL","C-TAL","dC-TAL","MD"))
KBR.sub.tal <- RenameIdents(object = KBR.sub.tal, 
                            "M-TAL" = "TAL",
                            "dM-TAL" = "dTAL",
                            "C-TAL" = "TAL",
                            "dC-TAL" = "dTAL",
                            "MD" = "TAL")
celltype <- levels(Idents(KBR.sub.tal))
cons.tal.DEGs <- getconservedDEGs(KBR.sub.tal, celltype)



###Conserved Degen State Markers - DCT
KBR.sub.dct <- subset(KBR.sub, idents = c("DCT1","DCT2","dDCT"))
KBR.sub.dct <- RenameIdents(object = KBR.sub.dct, 
                            "DCT1" = "DCT",
                            "DCT2" = "DCT")
KBR.sub.dct$subclass.l3 <- Idents(KBR.sub.dct)

celltype <- levels(Idents(KBR.sub.dct))
cons.dct.DEGs <- getconservedDEGs(KBR.sub.dct, celltype)




###Conserved Degen State Markers - CNT
KBR.sub.cnt <- subset(KBR.sub, idents = c("CNT","CNT-PC","dCNT"))
KBR.sub.cnt <- RenameIdents(object = KBR.sub.cnt, 
                            "CNT-PC" = "CNT")
celltype <- levels(Idents(KBR.sub.cnt))
cons.cnt.DEGs <- getconservedDEGs(KBR.sub.cnt, celltype)




###Conserved Degen State Markers - PC  (not including IMCD)
KBR.sub.pc <- subset(KBR.sub, idents = c("CCD-PC","OMCD-PC","dOMCD-PC"))
KBR.sub.pc <- RenameIdents(object = KBR.sub.pc, 
                           "CCD-PC" = "PC",
                           "OMCD-PC" = "PC",
                           "dOMCD-PC" = "dPC")
celltype <- levels(Idents(KBR.sub.pc))
cons.pc.DEGs <- getconservedDEGs(KBR.sub.pc, celltype)




###Conserved Degen State Markers - IC-A
KBR.sub.ic <- subset(KBR.sub, idents = c("CCD-IC-A","CNT-IC-A","dC-IC-A","OMCD-IC-A"))
KBR.sub.ic <- RenameIdents(object = KBR.sub.ic, 
                           "CCD-IC-A" = "IC-A",
                           "CNT-IC-A" = "IC-A",
                           "dC-IC-A" = "dIC-A",
                           "OMCD-IC-A" = "IC-A")
celltype <- levels(Idents(KBR.sub.ic))
cons.ic.DEGs <- getconservedDEGs(KBR.sub.ic, celltype)



###Conserved Degen State Markers - EC
KBR.sub.ec <- subset(KBR.sub, idents = c("EC-GC","EC-AEA","EC-DVR","EC-PTC","dEC-PTC",
                                         "EC-AVR","dEC"))
KBR.sub.ec <- RenameIdents(object = KBR.sub.ec, 
                           "EC-GC" = "EC",
                           "EC-AEA" = "EC",
                           "EC-DVR" = "EC",
                           "EC-PTC" = "EC",
                           "dEC-PTC" = "dEC",
                           "EC-AVR" = "EC",
                           "dEC" = "dEC")
celltype <- levels(Idents(KBR.sub.ec))
cons.ec.DEGs <- getconservedDEGs(KBR.sub.ec, celltype)



###Conserved Degen State Markers - VSMC/P
KBR.sub.smc <- subset(KBR.sub, idents = c("MC","REN","VSMC","VSMC/P","dVSMC"))
KBR.sub.smc <- RenameIdents(object = KBR.sub.smc, 
                            "MC" = "VSMC",
                            "REN" = "VSMC",
                            "VSMC" = "VSMC",
                            "VSMC/P" = "VSMC",
                            "dVSMC" = "dVSMC")
celltype <- levels(Idents(KBR.sub.smc))
cons.smc.DEGs <- getconservedDEGs(KBR.sub.smc, celltype)



###Conserved Degen State Markers - FIB - All
KBR.sub.fib <- subset(KBR.sub, idents = c("FIB","M-FIB","dM-FIB","dFIB"))
KBR.sub.fib <- RenameIdents(object = KBR.sub.fib, 
                            "M-FIB" = "FIB",
                            "dM-FIB" = "dFIB")
celltype <- levels(Idents(KBR.sub.fib))
cons.fib.DEGs <- getconservedDEGs(KBR.sub.fib, celltype)




###Combine markers for conserved gene set
DEG.list = list(cons.POD.DEGs,
                cons.pts13.DEGs,
                cons.dtl3.DEGs,
                cons.tal.DEGs,
                cons.dct.DEGs,
                cons.cnt.DEGs,
                cons.pc.DEGs,
                cons.ic.DEGs,
                cons.ec.DEGs,
                cons.smc.DEGs,
                cons.fib.DEGs)
celltype = 1:11

cons.DEGs <- do.call(rbind, lapply(celltype, function(ct) {
  print(paste("Running for cell type:", ct))
  
  cl.mark <- DEG.list[[ct]]
  cl.mark <- cl.mark[cl.mark$AKI_p_val < 0.05 & 
                       cl.mark$CKD_p_val < 0.05 & 
                       cl.mark$Ref_p_val < 0.05,]
  return(cl.mark)
  
}))




###Add in Degen score
degen.states <- c("dPOD", "dPT", "dDTL-ATL", "dTAL", "dDCT",  "dCNT", "dPC", "dIC-A",
                  "dEC", "dVSMC", "dFIB")
degen.gene.counts <- data.frame(cons.DEGs[cons.DEGs$cluster %in% degen.states,] %>%
                                  group_by(gene) %>%
                                  tally())
degen.gene.counts <- degen.gene.counts[order(-degen.gene.counts$n),]
degen.genes <- degen.gene.counts[degen.gene.counts$n >= 4, ]$gene #154
count.data = GetAssayData(KBR, slot = "data")
count.data = as.matrix(count.data[rownames(count.data) %in% toupper(degen.genes) ,])
combined_matrix = rbind(count.data,average=apply(count.data,2,mean))
cor_matrix = cor(t(combined_matrix))
cor_vector = cor_matrix[,dim(cor_matrix)[1]]
count.data_restricted = count.data[rownames(count.data) %in% names(cor_vector[cor_vector >= 0.1]),]
aaa = apply(count.data_restricted,2,mean)

KBR$degen.score <- aaa[rownames(KBR@meta.data)]



###Combine markers for expanded gene set
DEG.list = list(cons.POD.DEGs,
                cons.pts13.DEGs,
                cons.dtl3.DEGs,
                cons.tal.DEGs,
                cons.dct.DEGs,
                cons.cnt.DEGs,
                cons.pc.DEGs,
                cons.ic.DEGs,
                cons.ec.DEGs,
                cons.smc.DEGs,
                cons.fib.DEGs)
celltype = 1:11

cons.DEGs <- do.call(rbind, lapply(celltype, function(ct) {
  print(paste("Running for cell type:", ct))
  
  cl.mark <- DEG.list[[ct]]
  cl.mark <- subset(cl.mark, AKI_avg_log2FC >= 0.6 | CKD_avg_log2FC >= 0.6 | Ref_avg_log2FC >= 0.6)
  cl.mark <- subset(cl.mark, AKI_p_val_adj < 0.05 | CKD_p_val_adj < 0.05 | Ref_p_val_adj < 0.05)
  cl.mark <- cl.mark[cl.mark$max_pval < 0.05,]
  return(cl.mark)
  
}))




