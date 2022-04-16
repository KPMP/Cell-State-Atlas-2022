# scCv3 - Cycling State Gene Set and Score (Supplementary Table 15) --------------------------------------------------
library(Seurat)
library(SeuratDisk)
library(dplyr)
library(Matrix)
library(ggplot2)
library(viridis)
library(swne)

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


###Cycling Markers - PT
KBR.sub.pts13 <- subset(sc.KBR, idents = c("PT-S1/2","PT-S3","cycPT"))
KBR.sub.pts13 <- RenameIdents(object = KBR.sub.pts13, 
                              "PT-S1/2" = "PT",
                              "PT-S3" = "PT")

pts13.DEGs <- FindAllMarkers(KBR.sub.pts13, min.pct = 0.25, only.pos = TRUE, 
                             max.cells.per.ident = 300)




###Cycling Markers - EC
KBR.sub.ec <- subset(sc.KBR, idents = c("EC-GC","EC-AEA/DVR","EC-PTC",
                                        "EC-AVR","cycEC"))
KBR.sub.ec <- RenameIdents(object = KBR.sub.ec, 
                           "EC-GC" = "EC",
                           "EC-AEA/DVR" = "EC",
                           "EC-PTC" = "EC",
                           "EC-AVR" = "EC")
ec.DEGs <- FindAllMarkers(KBR.sub.ec, min.pct = 0.25, only.pos = TRUE, 
                          max.cells.per.ident = 300)




###Combine markers for conserved gene set
DEG.list = list(pts13.DEGs,
                ec.DEGs)
celltype = 1:2

cons.DEGs <- do.call(rbind, lapply(celltype, function(ct) {
  print(paste("Running for cell type:", ct))
  
  cl.mark <- DEG.list[[ct]]
  cl.mark <- cl.mark[cl.mark$p_val_adj < 0.05 & cl.mark$avg_log2FC > 0.6,]
  return(cl.mark)
  
}))

cyc.states <- c("cycPT","cycEC")
cyc.gene.counts <- data.frame(cons.DEGs[cons.DEGs$cluster %in% cyc.states,] %>%
                                group_by(gene) %>%
                                tally())
cyc.gene.counts <- cyc.gene.counts[order(-cyc.gene.counts$n),]





###Cycling state score
cyc.state.genes <- cyc.gene.counts[cyc.gene.counts$n == 2, ]$gene #101
count.data = GetAssayData(sc.KBR, slot = "data")
count.data = as.matrix(count.data[rownames(count.data) %in% toupper(cyc.state.genes) ,])
combined_matrix = rbind(count.data,average=apply(count.data,2,mean))
cor_matrix = cor(t(combined_matrix))
cor_vector = cor_matrix[,dim(cor_matrix)[1]]
count.data_restricted = count.data[rownames(count.data) %in% names(cor_vector[cor_vector >= 0.1]),]
aaa = apply(count.data_restricted,2,mean)
sc.KBR$cycling.score <- aaa[rownames(sc.KBR@meta.data)]







