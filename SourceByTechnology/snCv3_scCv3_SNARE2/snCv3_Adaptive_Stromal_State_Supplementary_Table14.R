# snCv3 - Adaptive Stromal State Gene Set and Score (Supplementary Table 14) -------------------------------------
library(Seurat)
library(SeuratDisk)
library(dplyr)
library(Matrix)
library(ggplot2)
library(viridis)
library(swne)

KBR <- LoadH5Seurat("Kidney_Healthy-Injury_Cell_Atlas_snCv3_Seurat_03282022.h5Seurat")
Idents(object = KBR) <- "clusters"
Idents(object = KBR) <- factor(Idents(object = KBR), levels = 1:119)
KBR.sub <- subset(KBR, idents = 1:100)



###Conserved Adaptive Stromal Markers
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
  "MAC-M2","cycMNP","MDC","cDC","pDC","ncMON","N","SC/NEU"))




###Conserved Markers of aFIB
KBR.sub.fib <- subset(KBR.sub, idents = c("FIB","aFIB"))
KBR.sub.fib$subclass.l3 <- Idents(KBR.sub.fib)
celltype <- levels(Idents(KBR.sub.fib))
cons.fib.as.DEGs <- getconservedDEGs(KBR.sub.fib, celltype)




###Conserved Markers of MYOF
KBR.sub.fib <- subset(KBR.sub, idents = c("MYOF","FIB"))
Idents(KBR.sub.fib) <- factor(Idents(KBR.sub.fib), levels = c("FIB", "MYOF"))
KBR.sub.fib$subclass.l3 <- Idents(KBR.sub.fib)
celltype <- levels(Idents(KBR.sub.fib))
cons.myof.as.DEGs <- getconservedDEGs(KBR.sub.fib, celltype)





###Conserved Markers for aFIB and MYOF combined
KBR.sub.fib <- subset(KBR.sub, idents = c("MYOF","FIB","aFIB"))
KBR.sub.fib <- RenameIdents(KBR.sub.fib, 'MYOF' = "aFIB")
Idents(KBR.sub.fib) <- factor(Idents(KBR.sub.fib), levels = c("FIB", "aFIB"))
KBR.sub.fib$subclass.l3 <- Idents(KBR.sub.fib)
celltype <- levels(Idents(KBR.sub.fib))
cons.aStr.DEGs <- getconservedDEGs(KBR.sub.fib, celltype)


###Conserved markers across aStr
###Combine markers
DEG.list = list(cons.aStr.DEGs)
celltype = 1

cons.DEGs <- do.call(rbind, lapply(celltype, function(ct) {
  print(paste("Running for cell type:", ct))
  
  cl.mark <- DEG.list[[ct]]
  cl.mark <- cl.mark[cl.mark$AKI_p_val < 0.05 & 
                       cl.mark$CKD_p_val < 0.05 & 
                       cl.mark$Ref_p_val < 0.05,]
  
  return(cl.mark)
  
}))



###Select conserved markers that are more specific to the altered state 
KBR.sub <- RenameIdents(KBR.sub, "MYOF" = "aFIB")
aFIB.markers <- FindMarkers(KBR.sub, ident.1 = "aFIB", features = aStr.genes,
                            only.pos = TRUE)




###Altered state score
alt.state.genes <- rownames(aFIB.markers[aFIB.markers$avg_log2FC > 0.6,]) #42 genes
count.data = GetAssayData(KBR, slot = "data")
count.data = as.matrix(count.data[rownames(count.data) %in% toupper(alt.state.genes) ,])
combined_matrix = rbind(count.data,average=apply(count.data,2,mean))
cor_matrix = cor(t(combined_matrix))
cor_vector = cor_matrix[,dim(cor_matrix)[1]]
count.data_restricted = count.data[rownames(count.data) %in% names(cor_vector[cor_vector >= 0.1]),]
aaa = apply(count.data_restricted,2,mean)
KBR$aStr.state.score <- aaa[rownames(KBR@meta.data)]




###Expanded adaptive -str state gene set
DEG.list = list(cons.fib.as.DEGs, cons.myof.as.DEGs)
celltype = 1:2

cons.DEGs <- do.call(rbind, lapply(celltype, function(ct) {
  print(paste("Running for cell type:", ct))
  
  cl.mark <- DEG.list[[ct]]
  cl.mark <- subset(cl.mark, AKI_avg_log2FC >= 0.6 | CKD_avg_log2FC >= 0.6 | Ref_avg_log2FC >= 0.6)
  cl.mark <- subset(cl.mark, AKI_p_val_adj < 0.05 | CKD_p_val_adj < 0.05 | Ref_p_val_adj < 0.05)
  cl.mark <- cl.mark[cl.mark$max_pval < 0.05,]
  return(cl.mark)
  
}))








