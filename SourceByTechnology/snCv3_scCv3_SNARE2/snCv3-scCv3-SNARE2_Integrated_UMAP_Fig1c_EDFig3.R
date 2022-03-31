# snCv3/scCv3/SNARE2 - Integrated UMAP (Fig 1c, ED Fig 3) ------------------------------
library(Signac)
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
set.seed(1234)

load("color_factors.robj")
load("sc-sn_int_color_factors.robj")
full.sc.l3.cols <- c(sc.l3.cols, hsc.l3.cols)



###Combined 10X/SNARE/sc10X embedding (Fig. 1c)
load ("Kidney_Healthy-Injury_Cell_Atlas_SNARE2-RNA-AC_Seurat_03282022.rda")
KSAC[['ATAC']] <- NULL
KSAC[['chromvar']] <- NULL
KSAC[['Activity']] <- NULL

snCv3 <- LoadH5Seurat("Kidney_Healthy-Injury_Cell_Atlas_snCv3_Seurat_03282022.h5Seurat")
scCv3 <- LoadH5Seurat("Kidney_Healthy-Injury_Cell_Atlas_scCv3_Seurat_03282022.h5Seurat")

combined <- merge(KSAC, y = c(snCv3,scCv3))

umap.coordinates <- rbind(Embeddings(object = KSAC, reduction = "umap.rna"),
                          Embeddings(object = snCv3, reduction = "umap"),
                          Embeddings(object = scCv3, reduction = "ref.umap"))
combined[["umap"]] <- CreateDimReducObject(embeddings = umap.coordinates, key = "umap_", assay = DefaultAssay(combined))

Idents(combined) <- "subclass.l3"
DimPlot(combined, reduction = "umap", label = FALSE, raster=FALSE,
        label.size = 4, repel = TRUE, shuffle = TRUE) + ggtitle("Subclass.l3"
        ) + scale_color_manual(values = alpha(full.sc.l3.cols[levels(Idents(combined))], 0.4), name = "Subclass.l3" 
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2))) + NoLegend()

Idents(combined) <- "assay"
combined <- RenameIdents(combined, '10X snRNA-Seq' = "10X snRNA-seq")
assay.cols <- c("#507EB3","gray","gray")
names(assay.cols) <- c("10X snRNA-seq","SNARE-Seq2","10X scRNA-seq")
DimPlot(combined, reduction = "umap", label = FALSE, raster=FALSE,
        label.size = 4, repel = TRUE, shuffle = TRUE) + ggtitle("Subclass.l3"
        ) + scale_color_manual(values = alpha(assay.cols[levels(Idents(combined))], 0.4), name = "Subclass.l3" 
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2))) + NoLegend()
assay.cols <- c("gray","#6F3980","gray")
names(assay.cols) <- c("10X snRNA-seq","SNARE-Seq2","10X scRNA-seq")
DimPlot(combined, reduction = "umap", label = FALSE, raster=FALSE,
        label.size = 4, repel = TRUE, shuffle = TRUE) + ggtitle("Subclass.l3"
        ) + scale_color_manual(values = alpha(assay.cols[levels(Idents(combined))], 0.4), name = "Subclass.l3" 
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2))) + NoLegend()
assay.cols <- c("gray","gray","#4D7880")
names(assay.cols) <- c("10X snRNA-seq","SNARE-Seq2","10X scRNA-seq")
DimPlot(combined, reduction = "umap", label = FALSE, raster=FALSE,
        label.size = 4, repel = TRUE, shuffle = TRUE) + ggtitle("Subclass.l3"
        ) + scale_color_manual(values = alpha(assay.cols[levels(Idents(combined))], 0.4), name = "Subclass.l3" 
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2))) + NoLegend()





###UMAP of Subclass.l3 (Extended Data Fig. 3a)
##Merge to concensus labels 
Idents(combined) <- "subclass.l3"
combined <- RenameIdents(combined, 
                         'PT-S1' = 'PT-S1/2',
                         'PT-S2' = 'PT-S1/2',
                         'NKT' = 'NKC/T',
                         'dPT/DTL' = 'dPT')
combined$subclass.l3 <- Idents(combined)
Idents(combined) <- "assay"
combined <- RenameIdents(combined, '10X snRNA-Seq' = "10X snRNA-seq")
combined$assay <- Idents(combined)

Idents(combined) <- "subclass.l3"
col.order <- c(
  "POD","dPOD","PEC","PT-S1/2","PT-S3","aPT","cycPT","dPT",
  "dPT/DTL","DTL2","DTL1","DTL3","dDTL3","ATL","dATL","aTAL1",
  "aTAL2","M-TAL","dM-TAL","C-TAL","dC-TAL","MD","DCT1","DCT2",
  "dDCT","cycDCT","CNT","CNT-PC","dCNT","cycCNT","CCD-PC","OMCD-PC",
  "dOMCD-PC","IMCD","dIMCD","CCD-IC-A","CNT-IC-A","dC-IC-A","OMCD-IC-A",
  "tPC-IC","IC-B","PapE","EC-GC","EC-AEA","EC-DVR","EC-PTC","dEC-PTC",
  "EC-AVR","dEC","cycEC","EC-LYM","MC","REN","VSMC","VSMC/P","dVSMC",
  "MYOF","cycMYOF","FIB","M-FIB","dM-FIB","aFIB","dFIB","B","PL","T",
  "NKC/T","cycNKC/T","MAST","MAC-M2","cycMNP","MDC","cDC","pDC","ncMON","N","SC/NEU")
Idents(combined) <- factor(Idents(combined), levels = col.order)
table(Idents(combined))

DimPlot(combined, reduction = "umap", label = TRUE, raster=FALSE,
        label.size = 4, repel = TRUE) + ggtitle("Subclass.l3"
        ) + scale_color_manual(values = alpha(full.sc.l3.cols[levels(Idents(combined))], 0.4), name = "Subclass.l3" 
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2))) + NoLegend()





###Plot QC stats per integrated subclass (Extended Data Fig. 3e)
col.order <- c(
  "POD","dPOD","PEC","PT-S1/2","PT-S3","aPT","cycPT","dPT",
  "dPT/DTL","DTL2","DTL1","DTL3","dDTL3","ATL","dATL","aTAL1",
  "aTAL2","M-TAL","dM-TAL","C-TAL","dC-TAL","MD","DCT1","DCT2",
  "dDCT","cycDCT","CNT","CNT-PC","dCNT","cycCNT","CCD-PC","OMCD-PC",
  "dOMCD-PC","IMCD","dIMCD","CCD-IC-A","CNT-IC-A","dC-IC-A","OMCD-IC-A",
  "tPC-IC","IC-B","PapE","EC-GC","EC-AEA","EC-DVR","EC-PTC","dEC-PTC",
  "EC-AVR","dEC","cycEC","EC-LYM","MC","REN","VSMC","VSMC/P","dVSMC",
  "MYOF","cycMYOF","FIB","M-FIB","dM-FIB","aFIB","dFIB","B","PL","T",
  "NKC/T","cycNKC/T","MAST","MAC-M2","cycMNP","MDC","cDC","pDC","ncMON","N","SC/NEU")

col.order <- col.order[col.order %in% levels(Idents(combined))]

load("sn-snare_int_color_factors.robj")
snare.patient.cols <- int.patient.cols
load("color_factors.robj")
load("sc-sn_int_color_factors.robj")

patient.cols <- c(int.patient.cols, snare.patient.cols)


Idents(combined) <- "assay" 
assay.cols <- c("#507EB3","#6F3980","#4D7880")
names(assay.cols) <- levels(Idents(combined))
Idents(combined) <- "subclass.l3" 

pdf(file='Subclass_Contributions_Barplot.pdf',width=8,height=20)
layout(matrix(c(1:8), nrow = 8, ncol = 1, byrow = TRUE))
prop1 <- prop.table(table(combined$assay, combined$subclass.l3), margin = 2)[,col.order]
barplot(prop1,main = "Assay Proportions",  names.arg = colnames(prop1), 
        cex.names = 0.5, las=2, col = as.character(assay.cols[rownames(prop1)]))
prop2 <- prop.table(table(combined$sex, combined$subclass.l3), margin = 2)[,col.order]
barplot(prop2, main = "Sex Proportions",  names.arg = colnames(prop2), 
        cex.names = 0.5, las=2, col = as.character(sex.cols[rownames(prop2)]))
prop3 <- prop.table(table(combined$condition.l1, combined$subclass.l3), margin = 2)[,col.order]
barplot(prop3, main = "Condition Proportions",  names.arg = colnames(prop3), 
        cex.names = 0.5, las=2, col = as.character(cond.l1.cols[rownames(prop3)]))
prop4 <- prop.table(table(combined$region.l1, combined$subclass.l3)[c("Cortex","Medulla","Papilla"),], margin = 2)[,col.order]
barplot(prop4, main = "Region Proportions",  names.arg = colnames(prop4), 
        cex.names = 0.5, las=2, col = as.character(region.l1.cols[rownames(prop4)]))
prop5 <- prop.table(table(combined$tissue_type, combined$subclass.l3), margin = 2)[,col.order]
barplot(prop5, main = "Tissue Type Proportions", names.arg = colnames(prop5),
        cex.names = 0.5, las=2, col = c("#E6ACAC","#805139","#4D8077"))
#Biopsy, Deceased Donor, Nephrectomy
prop6 <- prop.table(table(combined$patient, combined$subclass.l3), margin = 2)[,col.order]
barplot(prop6, main = "Patient Proportions", names.arg = colnames(prop6),
        cex.names = 0.5, las=2, col = as.character(patient.cols[rownames(prop6)]))
batch.entropy<-table(combined$patient, combined$subclass.l3)
entropy <- function(x) {
  value = entropy::entropy(x ,method='MM',unit='log2')/log2(length(x))
  return(value)
}
batch.entropy<-apply(batch.entropy,2,entropy)
barplot(batch.entropy[col.order], col = "gray", main = "Patient Entropy",  
        names.arg = col.order,las=2,
        cex.names = 0.5, ylim=c(0,1))
tab1 <- table(combined$subclass.l3)[col.order]
barplot(tab1, col = "gray", main = "Cluster Size", 
        names.arg = col.order,las=2, cex.names = 0.5)
dev.off()



###Correlation of Subclass.l3 between assays
sn.var <- VariableFeatures(snCv3)
VariableFeatures(combined) <- sn.var

snCv3 <- rownames(combined@meta.data[combined@meta.data$assay == "10X snRNA-seq",])
scCv3 <- rownames(combined@meta.data[combined@meta.data$assay == "10X scRNA-seq",])
SNARE2 <- rownames(combined@meta.data[combined@meta.data$assay == "SNARE-Seq2",])

snCv3 <- subset(combined, cells = snCv3)
snCv3 <- ScaleData(snCv3, features = sn.var)
scCv3 <- subset(combined, cells = scCv3)
scCv3 <- ScaleData(scCv3, features = sn.var)
SNARE2 <- subset(combined, cells = SNARE2)
SNARE2 <- ScaleData(SNARE2, features = sn.var)

##Use variable genes for correlation
ave.snCv3 <- AverageExpression(snCv3, features = sn.var, assays = "RNA",
                               slot = "scale.data")
ave.scCv3 <- AverageExpression(scCv3, features = sn.var, assays = "RNA",
                               slot = "scale.data")
ave.SNARE2 <- AverageExpression(SNARE2, features = sn.var, assays = "RNA",
                                slot = "scale.data")

ave.snsc.cor<-cor(ave.snCv3$RNA,ave.scCv3$RNA)
common.sn.sc <- intersect(rownames(ave.snsc.cor), colnames(ave.snsc.cor))
ave.snsc.cor<-ave.snsc.cor[common.sn.sc, common.sn.sc]

snsc.cor <- diag(ave.snsc.cor); names(snsc.cor) <- rownames(ave.snsc.cor)
snsc.cor[colnames(snsc.cor)==rownames(snsc.cor)]
ave.snsnare.cor<-cor(ave.snCv3$RNA,ave.SNARE2$RNA)
ave.snsnare.cor<-ave.snsnare.cor[rownames(ave.snsnare.cor) %in% colnames(ave.snsnare.cor), ]
snsnare.cor <- diag(ave.snsnare.cor)

merged <- merge(t(snsc.cor), t(snsnare.cor), all = TRUE)
rownames(merged) <- c("scCv3","SNARE2")
merged <- as.data.frame(t(merged))
merged <- merged[col.order,]
rownames(merged) <- col.order

library(gplots)
col2 <- colorRampPalette(c("#053061","#2166AC","#4393C3","#92C5DE","#D1E5F0","#FFFFFF","#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F"))
heatmap.2(t(merged),
          main = "Corr Values",
          density.info="none",
          trace="none",
          margins =c(8,4), 
          col=col2,
          dendrogram = "none",
          na.color="gray",
          Rowv = NULL,
          Colv = NULL,
          symm=F,symkey=F,symbreaks=T,
          #scale = "row",
          cexRow = 0.5,
          cexCol = 0.5 )       






###scCv3 UMAP Plots (Ext Data Fig 3b)
#Condition l2 
Idents(object = scCv3) <- "condition.l2"
DimPlot(scCv3, reduction = "ref.umap", label = FALSE, raster=FALSE,
        label.size = 4, repel = TRUE) + ggtitle("Condition"
        ) + scale_color_manual(values = alpha(cond.l2.cols[levels(Idents(scCv3))], 0.4), name = "Condition" 
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2))) 


#Sex 
Idents(object = scCv3) <- "sex"
DimPlot(scCv3, reduction = "ref.umap", label = FALSE, raster=FALSE,
        label.size = 4, repel = TRUE) + ggtitle("Sex"
        ) + scale_color_manual(values = alpha(sex.cols[levels(Idents(scCv3))], 0.4), name = "Sex" 
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2))) 


#Patient 
Idents(object = scCv3) <- "patient"
DimPlot(scCv3, reduction = "ref.umap", label = FALSE, raster=FALSE,
        label.size = 4, repel = TRUE) + ggtitle("Patient"
        ) + scale_color_manual(values = alpha(int.patient.cols[levels(Idents(scCv3))], 0.4), name = "Patient" 
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2))) 





###scCv3 Violin/Bar plots (Extended Data Fig 3d)
Idents(scCv3) <- "condition.l2"
Idents(scCv3) <- factor(Idents(scCv3), levels = c("Ref","AKI","DKD","H-CKD"))
scCv3$condition.l2 <- Idents(scCv3)
Idents(scCv3) <- "patient"

#Order by condition
order <- c(
  "1153EO1","1153EO2","1153EO3","1157EO1","1157EO2","1157EO3","1158EO1",
  "1158EO2","1158EO3","1162EO1","1162EO2","PRE018","PRE019","PRE027",
  "PRE038","PRE055","PRE062","PRE098","30-10018",
  "30-10034","30-10123","30-10125","32-10003","32-10034","32-10074","33-10005",
  "33-10006","34-10050","34-10184","34-10187","27-10039","28-10051","29-10006",
  "29-10010","29-10011","29-10012","29-10013","29-10016","31-10001","31-10035",
  "31-10040","31-10042","29-10008","31-10000","31-10013")
Idents(scCv3) <- factor(Idents(scCv3), levels = order)

VlnPlot(object = scCv3, features = c("nCount_RNA", "nFeature_RNA","percent.er","percent.mt"), 
        ncol = 1, pt.size = -1, cols = as.character(int.patient.cols[levels(Idents(scCv3))]))

pdf(file='scCv3_patient_Stats_Barplots.pdf',width=8,height=8)
layout(matrix(c(1,2,3), nrow = 3, ncol = 1, byrow = TRUE))
col.order <- order
row.order <- c("Ref","AKI","DKD","H-CKD")
prop1 <- prop.table(table(scCv3$condition.l2, scCv3$patient), margin = 2)[row.order,col.order]
barplot(prop1,main = "Condition Proportions", cex.names = 0.5, 
        names.arg = colnames(prop1),las=2, col = as.character(cond.l2.cols[rownames(prop1)]))
row.order <- c("POD","PEC","PT","DTL","ATL","TAL","DCT","CNT","PC","IC","PapE","EC","VSM/P",
               "FIB","IMM","NEU")
prop3 <- prop.table(table(scCv3$subclass.l1, scCv3$patient), margin = 2)[row.order,col.order]
barplot(prop3,main = "Subclass.l1 Proportions", cex.names = 0.5, 
        names.arg = colnames(prop3),las=2, col = as.character(sc.l1.cols[rownames(prop3)]))
tab1 <- table(scCv3$patient)[col.order]
barplot(tab1, col = "gray", main = "Cluster Size", cex.names = 0.5,las=2)
dev.off()






###SNARE2 UMAP Plots (ED Fig 3b)
load("sn-snare_int_color_factors.robj")

#Sex
Idents(object = KSAC) <- "sex"
DimPlot(KSAC, reduction = "umap.rna", label = FALSE, raster=FALSE, 
        label.size = 4, repel = TRUE) + ggtitle("Sex"
        ) + scale_color_manual(values = alpha(sex.cols[levels(Idents(KSAC))], 0.4), name = "Sex" 
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2))) 

#Patient
Idents(object = KSAC) <- "patient"
DimPlot(KSAC, reduction = "umap.rna", label = FALSE, raster=FALSE, 
        label.size = 4, repel = TRUE) + ggtitle("Patient"
        ) + scale_color_manual(values = alpha(int.patient.cols[levels(Idents(KSAC))], 0.4), name = "Patient" 
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2))) 

#Region.l2
Idents(object = KSAC) <- "region.l2"
DimPlot(KSAC, reduction = "umap.rna", label = FALSE, raster=FALSE, 
        label.size = 4, repel = TRUE) + ggtitle("Region.l2"
        ) + scale_color_manual(values = alpha(region.l2.cols[levels(Idents(object = KSAC))], 0.4), name = "Region" 
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2))) 






###SNARE2 QC stats per patient (ED Fig 3d)
load ("Kidney_Healthy-Injury_Cell_Atlas_SNARE2-RNA-AC_Seurat_03282022.rda")
Idents(KSAC) <- "patient"
col.order <- levels(Idents(KSAC))
pdf(file='SNARE2_Patient_Stats_Barplot.pdf',width=4,height=9)
layout(matrix(c(1,2,3), nrow = 3, ncol = 1, byrow = TRUE))
prop1 <- prop.table(table(KSAC$region.l2, KSAC$patient), margin = 2)[,col.order]
barplot(prop1,main = "Region Proportions",  names.arg = colnames(prop1), 
        cex.names = 0.5, las=2, col = as.character(region.l2.cols[rownames(prop1)]))
row.order <- c("POD","PEC","PT","DTL","ATL","TAL","DCT","CNT","PC","IC","EC","VSM/P",
               "FIB","IMM")
prop3 <- prop.table(table(KSAC$subclass.l1, KSAC$patient), margin = 2)[row.order,col.order]
barplot(prop3,main = "Subclass.l1 Proportions", cex.names = 0.5, 
        names.arg = colnames(prop3),las=2, col = as.character(sc.l1.cols[rownames(prop3)]))
tab1 <- table(KSAC$patient)[col.order]
barplot(tab1, col = "gray", main = "Sample Size", cex.names = 0.5,las=2)
dev.off()

VlnPlot(object = KSAC, features = c("nCount_RNA", "nFeature_RNA",
                                    "nCount_ATAC","TSS.enrichment"), 
        ncol = 1, pt.size = -1, cols = int.patient.cols[levels(Idents(KSAC))])





###SNARE2 Joint RNA/ATAC embedding (ED Fig 3c)
DefaultAssay(KSAC) <- "ATAC"

# build a joint neighbor graph using both assays (pre-computed)
KSAC <- FindMultiModalNeighbors(
  object = KSAC,
  reduction.list = list("pca", "lsi"), 
  dims.list = list(1:50, 2:30),
  modality.weight.name = "RNA.weight",
  verbose = TRUE
)

# build a joint UMAP visualization (pre-computed)
KSAC <- RunUMAP(
  object = KSAC,
  nn.name = "weighted.nn",
  assay = "RNA",
  verbose = TRUE
)

Idents(object = KSAC) <- "subclass.l3"
DimPlot(KSAC, reduction = "umap", label = TRUE, raster=FALSE,
        label.size = 4, repel = TRUE) + ggtitle("Subclass.l3"
        ) + NoLegend() + scale_color_manual(values = alpha(sc.l3.cols[levels(Idents(KSAC))], 0.4), name = "Subclass"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))

###

