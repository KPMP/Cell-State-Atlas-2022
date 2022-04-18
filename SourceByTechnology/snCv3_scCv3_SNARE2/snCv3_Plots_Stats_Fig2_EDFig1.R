# snCv3 - Plots and Stats (Fig 2 and ED Fig 1) -----------------------------------------------
library(Seurat)
library(SeuratDisk)
library(dplyr)
library(Matrix)
library(ggplot2)

KBR <- LoadH5Seurat("Kidney_Healthy-Injury_Cell_Atlas_snCv3_Seurat_03282022.h5Seurat")
load("color_factors.robj")



###UMAP Plots
##Generate plots for all clusters and annotated clusters
Idents(object = KBR) <- "clusters"
Idents(object = KBR) <- factor(Idents(object = KBR), levels = 1:119)
KBR.sub <- subset(KBR, idents = 1:100)

#Patient (Ext Data Fig. 1c)
Idents(object = KBR.sub) <- "patient"
DimPlot(KBR.sub, reduction = "umap", label = FALSE, raster=FALSE,
        label.size = 4, repel = TRUE) + ggtitle("Patient"
        ) + scale_color_manual(values = alpha(patient.cols[levels(Idents(KBR.sub))], 0.4), name = "Patient"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))


#Condition level 2 (Ext Data Fig. 1c)
Idents(object = KBR.sub) <- "condition.l2"
DimPlot(KBR.sub, reduction = "umap", label = FALSE, raster=FALSE,
        label.size = 4, repel = TRUE) + ggtitle("Condition"
        ) + scale_color_manual(values = alpha(cond.l2.cols[levels(Idents(KBR.sub))], 0.4), name = "Condition"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))


#State level 2 (Figure 2b)
Idents(object = KBR.sub) <- "state.l2"
DimPlot(KBR.sub, reduction = "umap", label = FALSE, raster=FALSE,
        label.size = 4, repel = TRUE) + ggtitle("State.l2"
        ) + scale_color_manual(values = alpha(state.l2.cols[levels(Idents(KBR.sub))], 0.4), name = "State.l2"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))


#Clusters (Ext Data Fig. 1a)
Idents(object = KBR.sub) <- "clusters"
Idents(object = KBR.sub) <- factor(Idents(object = KBR.sub), levels = 1:100)
DimPlot(KBR.sub, reduction = "umap", label = TRUE, raster=FALSE,
        label.size = 4, repel = TRUE) + ggtitle("Cluster"
        ) + NoLegend() + scale_color_manual(values = alpha(cl.cols[levels(Idents(object = KBR.sub))], 0.4), name = "Cluster"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))


#Subclass level 3 (Fig 2b)
Idents(object = KBR.sub) <- "subclass.l3"
DimPlot(KBR.sub, reduction = "umap", label = TRUE, raster=FALSE,
        label.size = 4, repel = TRUE) + ggtitle("Subclass.l3"
        ) + NoLegend() + scale_color_manual(values = alpha(sc.l3.cols[levels(Idents(KBR.sub))], 0.4), name = "Subclass"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))



#Sex (Ext Data Fig. 1c)
Idents(object = KBR.sub) <- "sex"
DimPlot(KBR.sub, reduction = "umap", label = FALSE, raster=FALSE,
        label.size = 4, repel = TRUE) + ggtitle("Sex"
        ) + scale_color_manual(values = alpha(sex.cols[levels(Idents(KBR.sub))], 0.4), name = "Sex"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))



#Region level 2 (Fig 2b, Ext Data Fig. 1c)
Idents(object = KBR.sub) <- "region.l2"
DimPlot(KBR.sub, reduction = "umap", label = FALSE, raster=FALSE,
        label.size = 4, repel = TRUE) + ggtitle("Region.l2"
        ) + scale_color_manual(values = alpha(region.l2.cols[levels(Idents(KBR.sub))], 0.4), name = "Region"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))



#snDROP Predicted ID (Ext Data Fig. 2b)
Idents(object = KBR.sub) <- "snDrop_predicted_subclass"
DimPlot(KBR.sub, reduction = "umap", label = TRUE, raster=FALSE,
        label.size = 4, repel = FALSE) + ggtitle("snDROP Predicted ID"
        ) + NoLegend() + scale_color_manual(values = alpha(snDrop.pred.cols[levels(Idents(KBR.sub))], 0.4), name = "snDROP Predicted ID"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))

FeaturePlot(KBR.sub, reduction = "umap", features = "snDrop_prediction_score", raster=FALSE)





###By Sample/patient Plots (Ext Data Fig. 1d)
Idents(KBR) <- "patient"

#Order by condition
order <- c("18-142", "18-162", "18-312", "3490", "3499", "3504", "3535", "3593", "3613",
           "KRP446", "KRP460", "KRP461", "KRP462", 
           "30-10034", "32-10003", "32-10034", "32-2", "33-10005", "33-10006", 
           "COV-02-0002", "COV-02-0010", "COV-02-0094 ", "COV-02-0096 ", 
           "29-10006", "29-10010", "29-10012", "29-10013", "31-10001", "31-10035","3487",
           "29-10008", "31-10000", "31-10006", "31-10013",
           "3477", "3479")
Idents(KBR) <- factor(Idents(KBR), levels = order)

#Violin plots
VlnPlot(object = KBR, features = c("nCount_RNA", "nFeature_RNA","percent.er","percent.mt"), 
        ncol = 1, pt.size = -1, cols = as.character(patient.cols[levels(Idents(KBR))]))

#Barplots
pdf(file='Patient_Stats_Barplots.pdf',width=8,height=12)
layout(matrix(c(1,2,3,4), nrow = 4, ncol = 1, byrow = TRUE))
col.order <- order
row.order <- c("Ref","AKI","COV-AKI","DKD","H-CKD")
prop1 <- prop.table(table(KBR$condition.l2, KBR$patient), margin = 2)[row.order,col.order]
barplot(prop1,main = "Condition Proportions", cex.names = 0.5, 
        names.arg = colnames(prop1),las=2, col = as.character(cond.l2.cols[rownames(prop1)]))
row.order <- c("C","C-M","M","P")
prop2 <- prop.table(table(KBR$region.l2, KBR$patient), margin = 2)[row.order,col.order]
barplot(prop2, main = "Region Proportions", cex.names = 0.5, 
        names.arg = colnames(prop2),las=2, col = as.character(region.l2.cols[rownames(prop2)]))
row.order <- c("POD","PEC","PT","DTL","ATL","TAL","DCT","CNT","PC","IC","PapE","EC","VSM/P",
               "FIB","IMM","NEU")
prop3 <- prop.table(table(KBR$subclass.l1, KBR$patient)[row.order,], margin = 2)[,col.order]
barplot(prop3,main = "Subclass.l1 Proportions", cex.names = 0.5, 
        names.arg = colnames(prop3),las=2, col = as.character(sc.l1.cols[rownames(prop3)]))
tab1 <- table(KBR$patient)[col.order]
barplot(tab1, col = "gray", main = "Cluster Size", cex.names = 0.5,las=2)
dev.off()




###By Cluster Plots (Ext Data Fig. 1f)
Idents(object = KBR.sub) <- "clusters"
Idents(object = KBR.sub) <- factor(Idents(object = KBR.sub), levels = 1:100)

#Violin plots
VlnPlot(object = KBR.sub, features = c("nCount_RNA", "nFeature_RNA","percent.er","percent.mt"), 
        ncol = 1, pt.size = -1, cols = cl.cols[levels(Idents(KBR.sub))])

#Barplots 1
pdf(file='Cluster_Barplots.pdf',width=8,height=12)
layout(matrix(c(1,2,3,4), nrow = 4, ncol = 1, byrow = TRUE))
col.order <- 1:100
prop1 <- prop.table(table(KBR.sub$condition.l2, KBR.sub$clusters), margin = 2)[,col.order]
barplot(prop1, main = "Condition Proportions", names.arg = colnames(prop1), las=2,
        cex.names = 0.5, col = as.character(cond.l2.cols[rownames(prop1)]))
prop2 <- prop.table(table(KBR.sub$patient, KBR.sub$clusters), margin = 2)[,col.order]
barplot(prop2, main = "Patient Proportions", names.arg = colnames(prop2),
        cex.names = 0.5, las=2, col = as.character(patient.cols[rownames(prop2)]))
prop3 <- prop.table(table(KBR.sub$region.l2, KBR.sub$clusters), margin = 2)[,col.order]
barplot(prop3,main = "Region Proportions", names.arg = colnames(prop3), las=2,
        cex.names = 0.5, col = as.character(region.l2.cols[rownames(prop3)]))
tab1 <- table(KBR.sub$clusters)[col.order]
barplot(tab1, col = "gray", main = "Cluster Size", 
        names.arg = col.order,las=2, cex.names = 0.5)
dev.off()


##






