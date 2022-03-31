# scCv3 - Degenerative State Gene Set and Score (Supplementary Table 12) --------------------------------------------------
library(Seurat)
library(SeuratDisk)
library(dplyr)
library(Matrix)
library(ggplot2)
library(forcats)
library(cowplot)
library(patchwork)

###Prepare seurat objects
sc.KBR <- LoadH5Seurat("Kidney_Healthy-Injury_Cell_Atlas_scCv3_Seurat_03282022.h5Seurat")
Idents(sc.KBR) <- "subclass.l3"
Idents(sc.KBR) <- factor(Idents(sc.KBR), levels = c(
  "POD","dPOD","PEC","PT-S1/2","PT-S3","aPT","cycPT","dPT",
  "dPT/DTL","DTL2","DTL1","DTL3","dDTL3","ATL","dATL","aTAL1",
  "aTAL2","M-TAL","dM-TAL","C-TAL","dC-TAL","MD","DCT1","DCT2",
  "dDCT","cycDCT","CNT","CNT-PC","dCNT","cycCNT","CCD-PC","OMCD-PC",
  "dOMCD-PC","IMCD","dIMCD","CCD-IC-A","CNT-IC-A","dC-IC-A","OMCD-IC-A",
  "tPC-IC","IC-B","PapE","EC-GC","EC-AEA","EC-DVR","EC-PTC","dEC-PTC",
  "EC-AVR","dEC","cycEC","EC-LYM","MC","REN","VSMC","VSMC/P","dVSMC",
  "MYOF","cycMYOF","FIB","M-FIB","dM-FIB","aFIB","dFIB","B","PL","T",
  "NKC/T","cycNKC/T","MAST","MAC-M2","cycMNP","MDC","cDC","pDC","ncMON","N","SC/NEU"))


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


###Conserved Degenerative Markers - POD
sc.KBR.pod <- subset(sc.KBR, idents = c("POD","dPOD"))
celltype <- levels(Idents(sc.KBR.pod))
cons.POD.DEGs <- getconservedDEGs(sc.KBR.pod, celltype)



###Conserved Degen Markers - PT S1-S3
sc.KBR.pts13 <- subset(sc.KBR, idents = c("PT-S1/2","PT-S3",
                                          "dPT"))
sc.KBR.pts13 <- RenameIdents(object = sc.KBR.pts13, 
                             "PT-S1/2" = "PT", 
                             "PT-S3" = "PT"
)
celltype <- levels(Idents(sc.KBR.pts13))
cons.pts13.DEGs <- getconservedDEGs(sc.KBR.pts13, celltype)



###Conserved Degen Markers - TAL
sc.KBR.tal <- subset(sc.KBR, idents = c("M-TAL","dM-TAL","C-TAL","dC-TAL","MD"))
sc.KBR.tal <- RenameIdents(object = sc.KBR.tal, 
                           "M-TAL" = "TAL",
                           "dM-TAL" = "dTAL",
                           "C-TAL" = "TAL",
                           "dC-TAL" = "dTAL",
                           "MD" = "TAL")
celltype <- levels(Idents(sc.KBR.tal))
cons.tal.DEGs <- getconservedDEGs(sc.KBR.tal, celltype)



###Conserved Degen Markers - DCT
sc.KBR.dct <- subset(sc.KBR, idents = c("DCT1","DCT2","dDCT"))
sc.KBR.dct <- RenameIdents(object = sc.KBR.dct, 
                           "DCT1" = "DCT",
                           "DCT2" = "DCT")
celltype <- levels(Idents(sc.KBR.dct))
cons.dct.DEGs <- getconservedDEGs(sc.KBR.dct, celltype)



###Conserved Degen Markers - CNT
sc.KBR.cnt <- subset(sc.KBR, idents = c("CNT","CNT-PC","dCNT"))
sc.KBR.cnt <- RenameIdents(object = sc.KBR.cnt, 
                           "CNT-PC" = "CNT")
celltype <- levels(Idents(sc.KBR.cnt))
cons.cnt.DEGs <- getconservedDEGs(sc.KBR.cnt, celltype)



###Conserved Degen Markers - PC  (not including IMCD)
sc.KBR.pc <- subset(sc.KBR, idents = c("CCD-PC","OMCD-PC","dOMCD-PC"))
sc.KBR.pc <- RenameIdents(object = sc.KBR.pc, 
                          "CCD-PC" = "PC",
                          "OMCD-PC" = "PC",
                          "dOMCD-PC" = "dPC")
celltype <- levels(Idents(sc.KBR.pc))
cons.pc.DEGs <- getconservedDEGs(sc.KBR.pc, celltype)



###Conserved Degen Markers - IC-A
sc.KBR.ic <- subset(sc.KBR, idents = c("CCD-IC-A","CNT-IC-A","dC-IC-A","OMCD-IC-A"))
sc.KBR.ic <- RenameIdents(object = sc.KBR.ic, 
                          "CCD-IC-A" = "IC-A",
                          "CNT-IC-A" = "IC-A",
                          "dC-IC-A" = "dIC-A",
                          "OMCD-IC-A" = "IC-A")
celltype <- levels(Idents(sc.KBR.ic))
cons.ic.DEGs <- getconservedDEGs(sc.KBR.ic, celltype)



###Conserved Degen Markers - VSMC/P
sc.KBR.smc <- subset(sc.KBR, idents = c("MC","VSMC","VSMC/P","dVSMC"))
sc.KBR.smc <- RenameIdents(object = sc.KBR.smc, 
                           "MC" = "VSMC",
                           "VSMC" = "VSMC",
                           "VSMC/P" = "VSMC",
                           "dVSMC" = "dVSMC")
celltype <- levels(Idents(sc.KBR.smc))
cons.smc.DEGs <- getconservedDEGs(sc.KBR.smc, celltype)



###Conserved Degen Markers - FIB 
sc.KBR.fib <- subset(sc.KBR, idents = c("FIB","M-FIB","dFIB"))
sc.KBR.fib <- RenameIdents(object = sc.KBR.fib, 
                           "M-FIB" = "FIB")
celltype <- levels(Idents(sc.KBR.fib))
cons.fib.DEGs <- getconservedDEGs(sc.KBR.fib, celltype)



###Combine markers
DEG.list = list(cons.POD.DEGs,
                cons.pts13.DEGs,
                cons.tal.DEGs,
                cons.dct.DEGs,
                cons.cnt.DEGs,
                cons.pc.DEGs,
                cons.ic.DEGs,
                cons.smc.DEGs,
                cons.fib.DEGs)
celltype = 1:9

cons.DEGs <- do.call(rbind, lapply(celltype, function(ct) {
  print(paste("Running for cell type:", ct))
  
  cl.mark <- DEG.list[[ct]]
  cl.mark <- cl.mark[cl.mark$AKI_p_val < 0.05 & 
                       cl.mark$CKD_p_val < 0.05 & 
                       cl.mark$Ref_p_val < 0.05,]
  top10 <- cl.mark %>% 
    mutate(avg_log2FC = (AKI_avg_log2FC + CKD_avg_log2FC + Ref_avg_log2FC) /3) %>% 
    group_by(cluster) %>% 
    top_n(n = 10, 
          wt = avg_log2FC)
  return(cl.mark)
  
}))




###Add in degen score
degen.states <- c("dPOD", "dPT", "dTAL", "dDCT",  "dCNT", "dPC", "dIC-A",
                  "dVSMC", "dFIB")
degen.gene.counts <- data.frame(cons.DEGs[cons.DEGs$cluster %in% degen.states,] %>%
                                  group_by(gene) %>%
                                  tally())
degen.gene.counts <- degen.gene.counts[order(-degen.gene.counts$n),]
degen.mark <- degen.gene.counts[degen.gene.counts$n >= 4, ]$gene
count.data = GetAssayData(sc.KBR, slot = "data")
count.data = as.matrix(count.data[rownames(count.data) %in% toupper(degen.mark) ,])
combined_matrix = rbind(count.data,average=apply(count.data,2,mean))
cor_matrix = cor(t(combined_matrix))
cor_vector = cor_matrix[,dim(cor_matrix)[1]]
count.data_restricted = count.data[rownames(count.data) %in% names(cor_vector[cor_vector >= 0.1]),]
aaa = apply(count.data_restricted,2,mean)
sc.KBR$degen.score <- aaa[rownames(sc.KBR@meta.data)]




###Combine markers for expanded gene set
DEG.list = list(cons.POD.DEGs,
                cons.pts13.DEGs,
                cons.tal.DEGs,
                cons.dct.DEGs,
                cons.cnt.DEGs,
                cons.pc.DEGs,
                cons.ic.DEGs,
                cons.smc.DEGs,
                cons.fib.DEGs)
celltype = 1:9

cons.DEGs <- do.call(rbind, lapply(celltype, function(ct) {
  print(paste("Running for cell type:", ct))
  
  cl.mark <- DEG.list[[ct]]
  cl.mark <- subset(cl.mark, AKI_avg_log2FC >= 0.6 | CKD_avg_log2FC >= 0.6 | Ref_avg_log2FC >= 0.6)
  cl.mark <- subset(cl.mark, AKI_p_val_adj < 0.05 | CKD_p_val_adj < 0.05 | Ref_p_val_adj < 0.05)
  cl.mark <- cl.mark[cl.mark$max_pval < 0.05,]
  return(cl.mark)
  
}))



