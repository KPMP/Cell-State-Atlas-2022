# This script is used to generate trajectory object for PT and TAL cells
# Qiwen Hu - 2021

library(Seurat)
library(slingshot)
library(pagoda2)
source("./Trajectory/util.func.R")

##############################
# Generate TAL trajectory
##############################
tal.seurat <- readRDS("./Trajectory/obj/tal.seurat.rds")
tal.p2 <- readRDS("./Trajectory/obj/tal.p2.removed.injury.rds") - # pagoda2 object

## remove M-tal
tal.seurat <- subset(tal.seurat, idents = c(26:30, 36:37))
tal.seurat <- RunUMAP(object = tal.seurat, reduction = "pca", dims = c(1:50), n.neighbors = 30L, min.dist = 0.01)

## Generate Slingshot object
tal.embeddings <- tal.seurat@reductions$umap@cell.embeddings
tal.annot <- setNames(tal.seurat@meta.data$subclass.l2, rownames(tal.seurat@meta.data))
tal.clusters <- setNames(tal.seurat@meta.data$clusters, rownames(tal.seurat@meta.data))
cellannot <- setNames(tal.seurat@meta.data$subclass.l2, rownames(tal.seurat@meta.data))

tal.sds <- slingshot(tal.embeddings, clusterLabels=tal.annot, approx_points=200, end.clus="ATL/TAL")
saveRDS(tal.sds, "./Trajectory/obj/tal.sds.removed.injury.rds")

## find trajectory associated genes
pt <- slingPseudotime(tal.sds, na=F)

mat <- tal.p2$counts
curve1.pt <- pt[, 1]
curve1.pt <- curve1.pt[names(curve1.pt) %in% rownames(mat)]
mat <- mat[rownames(mat) %in% names(curve1.pt), ]

mat <- t(mat[names(curve1.pt)[order(curve1.pt)], ])
ptime.genes <- test.associated.genes(curve1.pt, mat, n.cores=30)
saveRDS(ptime.genes, "./Trajectory/obj/tal.gam.assogene.rds")

##################################################################
# Generate PT trajectory
##################################################################
pt.seurat <- get(load("./Trajectory/obj/Kidney_10X-RNA_SNARE2-RNA_11302020_Integrated_Seurat_F_PT-only.rda"))

pt.cellannot <- factor(setNames(pt.seurat@meta.data$subclass.l3, names(Idents(pt.seurat))), levels=unique(pt.seurat@meta.data$subclass.l3))
pt.meta <- pt.seurat@meta.data
count <- pt.seurat@assays$RNA@counts[, colnames(pt.seurat@assays$RNA@counts) %in% rownames(pt.meta)]
pt.embedding <- pt.seurat@reductions$umap@cell.embeddings[rownames(pt.seurat@reductions$umap@cell.embeddings) %in% rownames(pt.meta), ]

pt.clusters <- pt.seurat$clusters[names(pt.seurat$clusters) %in% rownames(pt.embedding)]
pt.cellannot <- pt.cellannot[names(pt.cellannot) %in% rownames(pt.embedding)]
sds <- slingshot(pt.embedding, clusterLabels = pt.cellannot, approx_points=200, start.clus="PT-S1", end.clus="PT-S3")
saveRDS(sds, "./Trajectory/obj/pt.sds.rds")

## pagoda2 object
pt.p2 <- readRDS("./Trajectory/obj/pt.10X.snare.normal.altered.p2.rds")

# slingshot object
pt <- slingPseudotime(sds, na=F)

mat <- pt.p2.$counts
curve1.lambda <- sds@curves$curve1$lambda
curve1.pt <- pt[, 1]
curve1.pt <- curve1.pt[names(curve1.pt) %in% rownames(mat)]
mat <- mat[rownames(mat) %in% names(curve1.pt), ]

mat <- t(mat[names(curve1.pt)[order(curve1.pt)], ])
ptime.genes <- test.associated.genes(curve1.pt, mat, n.cores=30)
saveRDS(ptime.genes, "./Trajectory/obj/pt.gam.assogene.rds")






