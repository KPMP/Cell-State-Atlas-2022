# snCv3/scCv3/SNARE2 - Integrated UMAP (Fig 1) ------------------------------
library(Signac)
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
set.seed(1234)
load("color_factors.robj")
load("sc-sn_int_color_factors.robj")
full.sc.l3.cols <- c(sc.l3.cols, hsc.l3.cols)

###Combined 10X/SNARE/sc10X embedding
refquery <- LoadH5Seurat("Kidney_Healthy-Injury_Cell_Atlas_Integrated-snCv3-SNARE-RNA_Seurat_08032021.h5Seurat")
obj1 <- refquery
refquery <- LoadH5Seurat("Kidney_Healthy-Injury_Cell_Atlas_Integrated-snCv3-scCv3_Seurat_07302021.h5Seurat")
Idents(refquery) <- "id"
sc10x <- subset(refquery, idents = "query")
rm(refquery)
combined <- merge(obj1, y = sc10x)

umap.coordinates <- rbind(Embeddings(object = obj1, reduction = "ref.umap"),
                          Embeddings(object = sc10x, reduction = "ref.umap"))
combined[["ref.umap"]] <- CreateDimReducObject(embeddings = umap.coordinates, key = "refumap_", assay = DefaultAssay(combined))

Idents(combined) <- "subclass.l3"
DimPlot(combined, reduction = "ref.umap", label = FALSE, raster=FALSE,
        label.size = 4, repel = TRUE, shuffle = TRUE) + ggtitle("Subclass.l3"
        ) + scale_color_manual(values = alpha(full.sc.l3.cols[levels(Idents(combined))], 0.4), name = "Subclass.l3" 
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2))) + NoLegend()

Idents(combined) <- "assay"
combined <- RenameIdents(combined, '10X snRNA-Seq' = "10X snRNA-seq")
assay.cols <- c("#507EB3","gray","gray")
names(assay.cols) <- c("10X snRNA-seq","SNARE-Seq2","10X scRNA-seq")
DimPlot(combined, reduction = "ref.umap", label = FALSE, raster=FALSE,
        label.size = 4, repel = TRUE, shuffle = TRUE) + ggtitle("Subclass.l3"
        ) + scale_color_manual(values = alpha(assay.cols[levels(Idents(combined))], 0.4), name = "Subclass.l3" 
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2))) + NoLegend()
assay.cols <- c("gray","#6F3980","gray")
names(assay.cols) <- c("10X snRNA-seq","SNARE-Seq2","10X scRNA-seq")
DimPlot(combined, reduction = "ref.umap", label = FALSE, raster=FALSE,
        label.size = 4, repel = TRUE, shuffle = TRUE) + ggtitle("Subclass.l3"
        ) + scale_color_manual(values = alpha(assay.cols[levels(Idents(combined))], 0.4), name = "Subclass.l3" 
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2))) + NoLegend()
assay.cols <- c("gray","gray","#4D7880")
names(assay.cols) <- c("10X snRNA-seq","SNARE-Seq2","10X scRNA-seq")
DimPlot(combined, reduction = "ref.umap", label = FALSE, raster=FALSE,
        label.size = 4, repel = TRUE, shuffle = TRUE) + ggtitle("Subclass.l3"
        ) + scale_color_manual(values = alpha(assay.cols[levels(Idents(combined))], 0.4), name = "Subclass.l3" 
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2))) + NoLegend()




