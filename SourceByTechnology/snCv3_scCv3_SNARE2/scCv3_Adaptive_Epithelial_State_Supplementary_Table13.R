# scCv3 - Adaptive Epithelial State Gene Set and Score (Supplementary Table 13) ---------------------------------------------------------------
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




###Conserved Markers of aPT
KBR.sub.pt <- subset(sc.KBR, idents = c("PT-S1/2","PT-S3","aPT"))
KBR.sub.pt <- RenameIdents(object = KBR.sub.pt, 
                           "PT-S1/2" = "PT",
                           "PT-S3" = "PT")
table(Idents(KBR.sub.pt))
celltype <- levels(Idents(KBR.sub.pt))
cons.pt.as.DEGs <- getconservedDEGs(KBR.sub.pt, celltype)



###Conserved Markers of aTAL
KBR.sub.tal <- subset(sc.KBR, idents = c("aTAL1","M-TAL","C-TAL"))
KBR.sub.tal <- RenameIdents(object = KBR.sub.tal, 
                            "aTAL1" = "aTAL", 
                            "M-TAL" = "TAL",
                            "C-TAL" = "TAL")
Idents(KBR.sub.tal) <- factor(Idents(KBR.sub.tal), levels = c("TAL","aTAL"))
celltype <- levels(Idents(KBR.sub.tal))
cons.tal.as.DEGs <- getconservedDEGs(KBR.sub.tal, celltype)




###Conserved aEpi markers
###Combine markers
DEG.list = list(cons.pt.as.DEGs,
                cons.tal.as.DEGs)
celltype = 1:2

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

alt.states <- c("aPT","aTAL")
alt.gene.counts <- data.frame(cons.DEGs[cons.DEGs$cluster %in% alt.states,] %>%
                                group_by(gene) %>%
                                tally())
alt.gene.counts <- alt.gene.counts[order(-alt.gene.counts$n),]
common.alt.genes <- alt.gene.counts[alt.gene.counts$n == 2, ]$gene



###Select conserved markers that are more specific to the aEpi state 
sc.KBR <- RenameIdents(sc.KBR, "aPT" = "aTub")
sc.KBR <- RenameIdents(sc.KBR, "aTAL1" = "aTub")
sc.KBR <- RenameIdents(sc.KBR, "aTAL2" = "aTub")
aTub.markers <- FindMarkers(sc.KBR, ident.1 = "aTub", features = common.alt.genes,
                            only.pos = TRUE)




###Adaptive state score
adap.state.genes <- rownames(aTub.markers) #63 genes
adap.state.genes <- adap.state.genes[!adap.state.genes %in% "MALAT1"]
count.data = GetAssayData(sc.KBR, slot = "data")
count.data = as.matrix(count.data[rownames(count.data) %in% toupper(adap.state.genes) ,])
combined_matrix = rbind(count.data,average=apply(count.data,2,mean))
cor_matrix = cor(t(combined_matrix))
cor_vector = cor_matrix[,dim(cor_matrix)[1]]
count.data_restricted = count.data[rownames(count.data) %in% names(cor_vector[cor_vector >= 0.1]),]
aaa = apply(count.data_restricted,2,mean)
sc.KBR$aEpi.score <- aaa[rownames(sc.KBR@meta.data)]




###Expanded adaptive -epi state gene set
###Conserved markers across aPT and aTAL
###Combine markers

DEG.list = list(cons.pt.as.DEGs,
                cons.tal.as.DEGs)
celltype = 1:2

cons.DEGs <- do.call(rbind, lapply(celltype, function(ct) {
  print(paste("Running for cell type:", ct))
  
  cl.mark <- DEG.list[[ct]]
  cl.mark <- subset(cl.mark, AKI_avg_log2FC >= 0.6 | CKD_avg_log2FC >= 0.6 | Ref_avg_log2FC >= 0.6)
  cl.mark <- subset(cl.mark, AKI_p_val_adj < 0.05 | CKD_p_val_adj < 0.05 | Ref_p_val_adj < 0.05)
  cl.mark <- cl.mark[cl.mark$max_pval < 0.05,]
  return(cl.mark)
  
}))




