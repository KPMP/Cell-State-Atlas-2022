# SNARE2 - RNA Plots (Ext. Data Fig 4) --------------------------------------------------
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)

refquery <- LoadH5Seurat("Kidney_Healthy-Injury_Cell_Atlas_Integrated-snCv3-SNARE-RNA_Seurat_08032021.h5Seurat")
load("color_factors.robj")
load("sn-snare_int_color_factors.robj")

snare.cells <- rownames(refquery@meta.data[refquery@meta.data$id == "query",])

###UMAP Plots (Ext Data Fig 4a)
#Dataset Modalities
Idents(object = refquery) <- "id"
id.cols<-c("#36B3B3","#CC5252")
names(id.cols)<-c("reference","query")

#Sex
Idents(object = refquery) <- "sex"
DimPlot(refquery, reduction = "ref.umap", label = FALSE, raster=FALSE, cells = snare.cells,
        label.size = 4, repel = TRUE) + ggtitle("Sex"
        ) + scale_color_manual(values = alpha(sex.cols[levels(Idents(refquery))], 0.4), name = "Sex" 
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2))) 

#Patient
Idents(object = refquery) <- "patient"
DimPlot(refquery, reduction = "ref.umap", label = FALSE, raster=FALSE, cells = snare.cells,
        label.size = 4, repel = TRUE) + ggtitle("Patient"
        ) + scale_color_manual(values = alpha(int.patient.cols[levels(Idents(refquery))], 0.4), name = "Patient" 
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2))) 

#Region.l2
Idents(object = refquery) <- "region.l2"
DimPlot(refquery, reduction = "ref.umap", label = FALSE, raster=FALSE, cells = snare.cells,
        label.size = 4, repel = TRUE) + ggtitle("Region.l2"
        ) + scale_color_manual(values = alpha(region.l2.cols[levels(Idents(object = refquery))], 0.4), name = "Region" 
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2))) 

#Clusters
Idents(object = refquery) <- "clusters"
DimPlot(refquery, reduction = "ref.umap", label = TRUE, raster=FALSE, cells = snare.cells,
        label.size = 4, repel = TRUE) + ggtitle("Clusters"
        ) + scale_color_manual(values = alpha(cl.cols[levels(Idents(refquery))], 0.4), name = "Clusters" 
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2))) + NoLegend()



###Correlation of RNA clusters between snCv3 and SNARE (Ext Data Fig 4b)
Idents(object = refquery) <- "id"
ref <- subset(refquery, idents = "reference")
query <- subset(refquery, idents = "query")

Idents(ref) <- "clusters"
Idents(ref) <- factor(Idents(ref), levels = 1:100)

Idents(query) <- "clusters"
Idents(query) <- factor(Idents(query), levels = 1:100)

common.idents <- intersect(levels(Idents(ref)), levels(Idents(query)))
ref <- subset(ref, idents = common.idents)
query <- subset(query, idents = common.idents)

select.markers <- VariableFeatures(refquery)

ave.ref<-AverageExpression(ref, features = select.markers, assays = "RNA",
                           slot = "scale.data")
ave.query<-AverageExpression(query, features = select.markers, assays = "RNA",
                             slot = "scale.data")

library("corrplot")
ave.cor<-cor(ave.ref$RNA,ave.query$RNA)
col2 <- colorRampPalette(c("#053061","#2166AC","#4393C3","#92C5DE","#D1E5F0","#FFFFFF","#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F"))
corrplot(ave.cor, col = col2(500),  method="color",order = "original", cl.lim=c(-1,1))



###SNARE2 Barplot/Violin Plots (Ext Data Fig 4d)
###QC Plots
load("Kidney_Healthy-Injury_Cell_Atlas_SNARE2-RNA-AC_Seurat_08042021.rda")
KSAC

#Plot experiment-level QC stats per cluster
Idents(object = KSAC) <- "patient"
layout(matrix(c(1,2,3,4), nrow = 4, ncol = 1, byrow = TRUE))
prop1 <- prop.table(table(KSAC$patient, KSAC$clusters), margin = 2)[,col.order]
barplot(prop1,main = "Patient Proportions",  names.arg = colnames(prop1), 
        cex.names = 0.5, las=2, col = as.character(int.patient.cols[rownames(prop1)]))
prop2 <- prop.table(table(KSAC$region.l2, KSAC$clusters), margin = 2)[,col.order]
barplot(prop2,main = "Region Proportions", names.arg = colnames(prop2), las=2,
        cex.names = 0.5, col = as.character(region.l2.cols[rownames(prop2)]))
batch.entropy<-table(KSAC$experiment.short, KSAC$clusters)
entropy <- function(x) {
  value = entropy::entropy(x ,method='MM',unit='log2')/log2(length(x))
  return(value)
}
batch.entropy<-apply(batch.entropy,2,entropy)
barplot(batch.entropy[col.order], col = "gray", main = "Experiment Entropy",  
        names.arg = col.order,las=2,
        cex.names = 0.5, ylim=c(0,1))
tab1 <- table(KSAC$clusters)[col.order]
barplot(tab1, col = "gray", main = "Cluster Size", 
        names.arg = col.order,las=2, cex.names = 0.5)

#Violin plots - stats
Idents(KSAC) <- "clusters"
col.order <- 1:100
Idents(KSAC) <- factor(Idents(KSAC), levels = col.order)
VlnPlot(object = KSAC, features = c("nCount_RNA", "nFeature_RNA",
                                    "nCount_ATAC","TSS.enrichment",
                                    "nucleosome_signal","blacklist_fraction"), 
        ncol = 1, pt.size = -1, cols = cl.cols[levels(Idents(KSAC))])


###SNARE2-AC Cluster UMAP Plot (Ext Data Fig 4c)
#clusters
Idents(object = KSAC) <- "clusters"
DimPlot(KSAC, reduction = "umap.atac", label = TRUE, raster=FALSE,
        label.size = 4, repel = TRUE) + ggtitle("Clusters"
        ) + NoLegend() + scale_color_manual(values = alpha(cl.cols[levels(Idents(KSAC))], 0.4), name = "Clusters"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))




