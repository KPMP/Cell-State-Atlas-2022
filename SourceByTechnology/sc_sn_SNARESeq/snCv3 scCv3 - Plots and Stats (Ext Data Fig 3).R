# snCv3 scCv3 - Plots and Stats (Ext Data Fig 3)   --------------------------------------------------
library(Seurat)
library(SeuratDisk)
library(dplyr)
library(Matrix)
library(ggplot2)


###Prepare seurat objects
refquery <- LoadH5Seurat("Kidney_Healthy-Injury_Cell_Atlas_Integrated-snCv3-scCv3_Seurat_07302021.h5Seurat")
load("color_factors.robj")
load("sc-sn_int_color_factors.robj")


###Subset Single Cell data
Idents(object = refquery) <- "id"
sc.KBR <- subset(refquery, idents = "query")
Idents(sc.KBR) <- "subclass.l3"
Idents(sc.KBR) <- factor(Idents(sc.KBR), levels = c(
  "POD","dPOD","PEC","PT-S1/2","PT-S3","aPT","cycPT","dPT",
  "dPT/DTL","DTL2","DTL1","DTL3","dDTL3","ATL","dATL","aTAL1",
  "aTAL2","M-TAL","dM-TAL","C-TAL","dC-TAL","MD","DCT1","DCT2",
  "dDCT","cycDCT","CNT","CNT-PC","dCNT","cycCNT","CCD-PC","OMCD-PC",
  "dOMCD-PC","IMCD","dIMCD","CCD-IC-A","CNT-IC-A","dC-IC-A","OMCD-IC-A",
  "tPC-IC","IC-B","PapE","EC-GC","EC-AEA/DVR","EC-PTC","dEC-PTC",
  "EC-AVR","dEC","cycEC","EC-LYM","MC","REN","VSMC","VSMC/P","dVSMC",
  "MYOF","cycMYOF","FIB","M-FIB","dM-FIB","aFIB","dFIB","B","PL","T",
  "NKC/T","cycNKC/T","MAST","MAC-M2","cycMNP","MDC","cDC","pDC","ncMON","N","SC/NEU"))

#Update meta
ER.GO <- read.delim("Cytosolic_Ribosome_Gene_Ontology.txt",sep="\t",header=TRUE,row.names=1)
ER.GO <- as.character(ER.GO$Symbol)
ER.GO <- intersect(ER.GO, rownames(sc.KBR))
sc.KBR[["percent.er"]] <- PercentageFeatureSet(sc.KBR, features = ER.GO)
sc.KBR[["percent.mt"]] <- PercentageFeatureSet(sc.KBR, pattern = "^MT-")




###UMAP Plots

#Dataset Modalities (Ext Data Fig 3b)
Idents(object = refquery) <- "id"
id.cols<-c("#36B3B3","#CC5252")
names(id.cols)<-c("reference","query")
DimPlot(refquery, reduction = "ref.umap", label = FALSE, raster=FALSE,
        label.size = 4, repel = TRUE) + ggtitle("Dataset"
        ) + scale_color_manual(values = alpha(id.cols[levels(Idents(refquery))], 0.4), name = "Dataset" 
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2))) #+ NoLegend()
DimPlot(refquery, reduction = "umap", label = FALSE, raster=FALSE,
        label.size = 4, repel = TRUE) + ggtitle("Dataset"
        ) + scale_color_manual(values = alpha(id.cols[levels(Idents(refquery))], 0.4), name = "Dataset" 
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2))) #+ NoLegend()


id.cols <- c("#507EB3","gray")
names(id.cols) <- c("reference","query")
DimPlot(refquery, reduction = "umap", label = FALSE, raster=FALSE,
        label.size = 4, repel = TRUE) + ggtitle("Dataset"
        ) + scale_color_manual(values = alpha(id.cols[levels(Idents(refquery))], 0.4), name = "Dataset" 
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2))) #+ NoLegend()


id.cols <- c("gray","#4D7880")
names(id.cols) <- c("reference","query")
DimPlot(refquery, reduction = "umap", label = FALSE, raster=FALSE,
        label.size = 4, repel = TRUE) + ggtitle("Dataset"
        ) + scale_color_manual(values = alpha(id.cols[levels(Idents(refquery))], 0.4), name = "Dataset" 
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2))) #+ NoLegend()



#subclass.l3 (Ext Data Fig 3b, d)
Idents(refquery) <- "subclass.l3"
DimPlot(refquery, reduction = "ref.umap", label = TRUE,  raster=FALSE,
        label.size = 4, repel = TRUE) + ggtitle("Subclass_l3"
        ) + scale_color_manual(values = alpha(hsc.l3.cols[levels(Idents(refquery))], 0.4), name = "Subclass_l3" 
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2))) + NoLegend()
DimPlot(refquery, reduction = "umap", label = TRUE,  raster=FALSE,
        label.size = 4, repel = TRUE) + ggtitle("Subclass_l3"
        ) + scale_color_manual(values = alpha(hsc.l3.cols[levels(Idents(refquery))], 0.4), name = "Subclass_l3" 
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2))) + NoLegend()

scells <- rownames(refquery@meta.data[refquery@meta.data$id == "query",])
DimPlot(refquery, reduction = "ref.umap", label = TRUE, cells = scells, raster=FALSE,
        label.size = 4, repel = TRUE) + ggtitle("Subclass_l3"
        ) + scale_color_manual(values = alpha(hsc.l3.cols[levels(Idents(refquery))], 0.4), name = "Subclass_l3" 
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2))) + NoLegend()



#Condition l2 (Ext Data Fig 3d)
Idents(object = refquery) <- "condition.l2"
DimPlot(refquery, reduction = "ref.umap", label = FALSE, cells = scells, raster=FALSE,
        label.size = 4, repel = TRUE) + ggtitle("Condition"
        ) + scale_color_manual(values = alpha(cond.l2.cols[levels(Idents(refquery))], 0.4), name = "Condition" 
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2))) 


#Sex (Ext Data Fig 3d)
Idents(object = refquery) <- "sex"
DimPlot(refquery, reduction = "ref.umap", label = FALSE, cells = scells, raster=FALSE,
        label.size = 4, repel = TRUE) + ggtitle("Sex"
        ) + scale_color_manual(values = alpha(sex.cols[levels(Idents(refquery))], 0.4), name = "Sex" 
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2))) 


#Patient (Ext Data Fig 3d)
Idents(object = refquery) <- "patient"
DimPlot(refquery, reduction = "ref.umap", label = FALSE, cells = scells, raster=FALSE,
        label.size = 4, repel = TRUE) + ggtitle("Patient"
        ) + scale_color_manual(values = alpha(int.patient.cols[levels(Idents(refquery))], 0.4), name = "Patient" 
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2))) 





###Plot Subclass level QC stats (Ext Data Fig 3e)
order <- c("POD","dPOD","PEC","PT-S1/2","PT-S3","aPT","cycPT","dPT",
           "dPT/DTL","DTL2","DTL1","DTL3","dDTL3","ATL","dATL","aTAL1",
           "aTAL2","M-TAL","dM-TAL","C-TAL","dC-TAL","MD","DCT1","DCT2",
           "dDCT","cycDCT","CNT","CNT-PC","dCNT","cycCNT","CCD-PC","OMCD-PC",
           "dOMCD-PC","IMCD","dIMCD","CCD-IC-A","CNT-IC-A","dC-IC-A","OMCD-IC-A",
           "tPC-IC","IC-B","PapE","EC-GC","EC-AEA/DVR","EC-PTC","dEC-PTC",
           "EC-AVR","dEC","cycEC","EC-LYM","MC","REN","VSMC","VSMC/P","dVSMC",
           "MYOF","cycMYOF","FIB","M-FIB","dM-FIB","aFIB","dFIB","B","PL","T",
           "NKC/T","cycNKC/T","MAST","MAC-M2","cycMNP","MDC","cDC","pDC","ncMON","N","SC/NEU")

Idents(object = sc.KBR) <- "subclass.l3"
Idents(object = sc.KBR) <- factor(Idents(object = sc.KBR), levels = order)

#Violin plots
VlnPlot(object = sc.KBR, features = c("nCount_RNA", "nFeature_RNA","percent.er","percent.mt","degen.score"), 
        ncol = 1, pt.size = -1, cols = hsc.l3.cols[levels(Idents(sc.KBR))])

#Barplots - 1
layout(matrix(c(1,2,3), nrow = 3, ncol = 1, byrow = TRUE))
col.order <- levels(Idents(sc.KBR))
prop2 <- prop.table(table(sc.KBR$patient, sc.KBR$subclass.l3), margin = 2)[,col.order]
barplot(prop2, main = "Patient Proportions", names.arg = colnames(prop2),
        cex.names = 0.5, las=2, col = as.character(int.patient.cols[rownames(prop2)]))
batch.entropy<-table(sc.KBR$experiment, sc.KBR$subclass.l3)
entropy <- function(x) {
  value = entropy::entropy(x ,method='MM',unit='log2')/log2(length(x))
  return(value)
}
batch.entropy<-apply(batch.entropy,2,entropy)
barplot(batch.entropy[col.order], col = "gray", main = "Experiment Entropy",  
        names.arg = col.order,las=2,
        cex.names = 0.5, ylim=c(0,1))
tab1 <- table(sc.KBR$subclass.l3)[col.order]
barplot(tab1, col = "gray", main = "Cluster Size", 
        names.arg = col.order,las=2, cex.names = 0.5)



#Barplots - 2
layout(matrix(c(1,2,3), nrow = 3, ncol = 1, byrow = TRUE))
col.order <- levels(Idents(sc.KBR))
prop1 <- prop.table(table(sc.KBR$condition.l2, sc.KBR$subclass.l3), margin = 2)[,col.order]
barplot(prop1, main = "Condition Proportions", names.arg = colnames(prop1), las=2,
        cex.names = 0.5, col = as.character(cond.l2.cols[rownames(prop1)]))
batch.entropy<-table(sc.KBR$condition.l2, sc.KBR$subclass.l3)
entropy <- function(x) {
  value = entropy::entropy(x ,method='MM',unit='log2')/log2(length(x))
  return(value)
}
batch.entropy<-apply(batch.entropy,2,entropy)
barplot(batch.entropy[col.order], col = "gray", main = "Condition Entropy",  
        names.arg = col.order,las=2,
        cex.names = 0.5, ylim=c(0,1))
tab1 <- table(sc.KBR$subclass.l3)[col.order]
barplot(tab1, col = "gray", main = "Cluster Size", 
        names.arg = col.order,las=2, cex.names = 0.5)


###Plot UMAP of Original scCv3 Annotations (Ext Data Fig 3a)
load("Premiere_Biopsy_Integrated_Seurat_03102021.rda") # not provided
sc.meta <- as.data.frame(readRDS("Premiere_Metadata_InjuryScore_New_032321.RDS"))
sc.meta$umap1 <- as.numeric(sc.meta$umap1)
sc.meta$umap2 <- as.numeric(sc.meta$umap2)
umap.coordinates <- as.matrix(sc.meta[,c("umap1","umap2")])
sc.KBR <- subset(sc.KBR, cells = rownames(umap.coordinates))
umap.coordinates <- umap.coordinates[colnames(sc.KBR),] 
sc.KBR[["sc.umap"]] <- CreateDimReducObject(embeddings = umap.coordinates, key = "scumap_", assay = DefaultAssay(sc.KBR))

#subclass.l2
Idents(sc.KBR) <- "sc.subclass.l2"
levels(Idents(sc.KBR))
new <- c("#B36236","#80784D","#654D80","#E68AC9","#AB67E6","#654D80","#B3369F","#E66A45","#36B366","#E6456A","#7F36B3","#ACE6D8",
         "#a6609b","#B37936","#39807A")
names(new) <- c("dPT/DTL","PT/DTL","T-CYT","MAC-M2","NK1","NK2","NKT","PT-S1/S2","PC","MAC-M1","T-REG","IC-A",
                "cycT","cycEPI","dIC-A")
hsc.l3.cols <- c(hsc.l3.cols, new)

DimPlot(sc.KBR, reduction = "raji.umap", label = TRUE,  raster=FALSE,
        label.size = 4, repel = TRUE) + ggtitle("Subclass_l2"
        ) + scale_color_manual(values = alpha(hsc.l3.cols[levels(Idents(sc.KBR))], 0.4), name = "Subclass_l2" 
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2))) + NoLegend()





###Correlation of subclasses between scCv3 and snCv3 (Ext Data Fig 3c)
Idents(object = refquery) <- "id"
ref <- subset(refquery, idents = "reference")
query <- subset(refquery, idents = "query")
order <- c("POD","dPOD","PEC","PT-S1/2","PT-S3","aPT","cycPT","dPT",
           "dPT/DTL","DTL2","DTL1","DTL3","dDTL3","ATL","dATL","aTAL1",
           "aTAL2","M-TAL","dM-TAL","C-TAL","dC-TAL","MD","DCT1","DCT2",
           "dDCT","cycDCT","CNT","CNT-PC","dCNT","cycCNT","CCD-PC","OMCD-PC",
           "dOMCD-PC","IMCD","dIMCD","CCD-IC-A","CNT-IC-A","dC-IC-A","OMCD-IC-A",
           "tPC-IC","IC-B","PapE","EC-GC","EC-AEA/DVR","EC-PTC","dEC-PTC",
           "EC-AVR","dEC","cycEC","EC-LYM","MC","REN","VSMC","VSMC/P","dVSMC",
           "MYOF","cycMYOF","FIB","M-FIB","dM-FIB","aFIB","dFIB","B","PL","T",
           "NKC/T","cycNKC/T","MAST","MAC-M2","cycMNP","MDC","cDC","pDC","ncMON","N","SC/NEU")

Idents(ref) <- "subclass.l3"
Idents(ref) <- factor(Idents(ref), levels = order)
Idents(query) <- "subclass.l3"
Idents(query) <- factor(Idents(query), levels = order)

#Correlation of subclasses that have >10 sc or sn
tab <- table(query$subclass.l3,query$id)
common.idents <- names(tab[tab[,1] > 10,])
query <- subset(query, idents = common.idents)

#Use variable genes for correlation
select.markers <- VariableFeatures(refquery)
ave.ref<-AverageExpression(ref, features = select.markers, assays = "RNA",
                           slot = "scale.data")
ave.query<-AverageExpression(query, features = select.markers, assays = "RNA",
                             slot = "scale.data")

library("corrplot")
ave.cor<-cor(cbind(ave.ref$RNA,ave.query$RNA))
ave.cor<-ave.cor[1:74,75:136]
col2 <- colorRampPalette(c("#053061","#2166AC","#4393C3","#92C5DE","#D1E5F0","#FFFFFF","#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F"))
corrplot(ave.cor, col = col2(500),  method="color",order = "original",
         hclust.method = "ward.D", cl.lim=c(-1,1))


###









