# snCv3/scCv3 - Marker Genes ----------------------------------------------
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)

####Load objects object
scCv3 <- LoadH5Seurat("Kidney_Healthy-Injury_Cell_Atlas_scCv3_Seurat_03282022.h5Seurat")
snCv3 <- LoadH5Seurat("Kidney_Healthy-Injury_Cell_Atlas_snCv3_Seurat_03282022.h5Seurat")

###Find Subclass.l3 markers
subclass.l3.markers.sn <- FindAllMarkers(snCv3, features = rownames(snCv3),
                                         only.pos = TRUE, max.cells.per.ident = 2000,
                                         logfc.threshold = 0.25, min.pct = 0.1, return.thresh = 0.05)

subclass.l3.markers.sc <- FindAllMarkers(scCv3, features = rownames(scCv3),
                                         only.pos = TRUE, max.cells.per.ident = 2000,
                                         logfc.threshold = 0.25, min.pct = 0.1, return.thresh = 0.05)

###

