# scCv3 - Adaptive Stromal State Gene Set and Score (Supplementary Table 14) ---------------------------------------------------------------
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


#Load Functions
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



###Conserved Markers of aFIB
sc.KBR.fib <- subset(sc.KBR, idents = c("FIB","aFIB"))
sc.KBR.fib$subclass.l3 <- Idents(sc.KBR.fib)
celltype <- levels(Idents(sc.KBR.fib))
cons.fib.as.DEGs <- getconservedDEGs(sc.KBR.fib, celltype)




###Markers of MYOF
sc.KBR.fib <- subset(sc.KBR, idents = c("FIB","MYOF"))
Idents(sc.KBR.fib) <- factor(Idents(sc.KBR.fib), levels = c("FIB","MYOF"))
cons.myof.as.DEGs <- FindAllMarkers(sc.KBR.fib, min.pct = 0.25, only.pos = TRUE, 
                                    max.cells.per.ident = 300)



###Conserved Markers for aFIB and MYOF combined
sc.KBR.fib <- subset(sc.KBR, idents = c("MYOF","FIB","aFIB"))
sc.KBR.fib <- RenameIdents(sc.KBR.fib, 'MYOF' = "aFIB")
Idents(sc.KBR.fib) <- factor(Idents(sc.KBR.fib), levels = c("FIB", "aFIB"))
sc.KBR.fib$subclass.l3 <- Idents(sc.KBR.fib)
celltype <- levels(Idents(sc.KBR.fib))
cons.aStr.DEGs <- getconservedDEGs(sc.KBR.fib, celltype)




###Conserved markers across aStr
###Combine markers
DEG.list = list(cons.aStr.DEGs)
celltype = 1

cons.DEGs <- do.call(rbind, lapply(celltype, function(ct) {
  print(paste("Running for cell type:", ct))
  
  cl.mark <- DEG.list[[ct]]
  cl.mark <- cl.mark[cl.mark$AKI_p_val < 0.05 & 
                       #cl.mark$CKD_p_val < 0.05 & 
                       cl.mark$Ref_p_val < 0.05,]
  
  return(cl.mark)
  
}))

alt.states <- c("aFIB")
alt.gene.counts <- data.frame(cons.DEGs[cons.DEGs$cluster %in% alt.states,] %>%
                                group_by(gene) %>%
                                tally())
alt.gene.counts <- alt.gene.counts[order(-alt.gene.counts$n),]
aStr.genes <- alt.gene.counts$gene



###Select conserved markers that are more specific to the altered state 
Idents(sc.KBR) <- "subclass.l3"
sc.KBR <- RenameIdents(sc.KBR, "MYOF" = "aFIB")
myof.as.DEGs <- subset(cons.myof.as.DEGs, avg_log2FC >= 0.6)
myof.as.DEGs <- subset(myof.as.DEGs, p_val_adj < 0.05)
aFIB.markers <- FindMarkers(sc.KBR, ident.1 = "aFIB", features = c(aStr.genes,myof.as.DEGs$gene),
                            only.pos = TRUE)




###Altered state score
alt.state.genes <- rownames(aFIB.markers[aFIB.markers$avg_log2FC > 0.6,]) #23 genes
count.data = GetAssayData(sc.KBR, slot = "data")
count.data = as.matrix(count.data[rownames(count.data) %in% toupper(alt.state.genes) ,])
combined_matrix = rbind(count.data,average=apply(count.data,2,mean))
cor_matrix = cor(t(combined_matrix))
cor_vector = cor_matrix[,dim(cor_matrix)[1]]
count.data_restricted = count.data[rownames(count.data) %in% names(cor_vector[cor_vector >= 0.1]),]
aaa = apply(count.data_restricted,2,mean)
sc.KBR$aStr.state.score <- aaa[rownames(sc.KBR@meta.data)]




###Expanded adaptive -epi state gene set
###Conserved markers across aPT and aTAL
###Combine markers

fib.as.DEGs <- subset(cons.fib.as.DEGs, AKI_avg_log2FC >= 0.6 | CKD_avg_log2FC >= 0.6 | Ref_avg_log2FC >= 0.6)
fib.as.DEGs <- subset(fib.as.DEGs, AKI_p_val < 0.01 | CKD_p_val < 0.01 | Ref_p_val < 0.01)
myof.as.DEGs <- subset(cons.myof.as.DEGs, avg_log2FC >= 0.6)
myof.as.DEGs <- subset(myof.as.DEGs, p_val_adj < 0.05)
alt.states <- c("aFIB", "MYOF")
cons.DEGs <- rbind(fib.as.DEGs[fib.as.DEGs$cluster %in% alt.states,c(19,18)],
                   myof.as.DEGs[myof.as.DEGs$cluster %in% alt.states,c(6,7)])
alt.gene.counts <- data.frame(cons.DEGs[cons.DEGs$cluster %in% alt.states,] %>%
                                group_by(gene) %>%
                                tally())
alt.gene.counts <- alt.gene.counts[order(-alt.gene.counts$n),]




