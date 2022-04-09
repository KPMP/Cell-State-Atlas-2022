# This script is used to generate trajectory object for PT and TAL cells
# Qiwen Hu - 2021

library(Seurat)
library(slingshot)
library(pagoda2)
library(WGCNA)
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

########################################
## module detection
########################################
tal.sds <- readRDS("./Trajectory/obj/tal.sds.removed.injury.rds")
  
tal.mat <- t(tal.p2$counts)
cellWeights <- slingCurveWeights(tal.sds)
curve1.lambda <- tal.sds@curves$curve1$lambda
curve1.tal <- tal.pt[, 1]
names(curve1.tal) <- rownames(tal.pt)
curve1.tal <- curve1.tal[names(curve1.tal) %in% colnames(tal.mat)]
curve1.tal <- curve1.tal[!is.na(curve1.tal)]

tal.sig.de.genes <- read.table("./Trajectory/obj/tal.sig.de.gene.txt", header=TRUE, sep="\t")
tal.mat <- tal.mat[rownames(tal.mat) %in% tal.sig.de.genes$gene, ]
tal.mat <- tal.mat[, colnames(tal.mat) %in% names(curve1.tal)]

### smoothing
gene_mat <- zoo::rollapply(as.matrix(t(tal.mat)), 30, mean, 
                           align='left',partial=T) %>% apply(2, rescale.and.center, max.quantile=1-1e-3) %>% t
colnames(gene_mat) <- colnames(tal.mat)
saveRDS(gene_mat, "/home/qiwenhu/hubmap/kidney/figs/tables/tal.sig.cond.gene_mat.1.rds")
powers <- c(c(1:10), seq(from = 12, to=20, by=2));
sft=pickSoftThreshold(t(gene_mat),
                      dataIsExpr = TRUE,powerVector = powers,
                      corFnc = cor, corOptions = list(use = 'p'),networkType = "signed")

softPower = 10
#calclute the adjacency matrix
adj= adjacency(t(gene_mat), type = "signed", power = softPower);
#turn adjacency matrix into a topological overlap matrix (TOM) to minimize the effects of noise and spurious associations
TOM=TOMsimilarityFromExpr(t(gene_mat),networkType = "signed", TOMType = "signed", power = softPower)
colnames(TOM) <- rownames(TOM) <- rownames(gene_mat)
dissTOM=1-TOM

geneTree = flashClust::flashClust(as.dist(dissTOM),method="average")
#plot the resulting clustering tree (dendrogram)
plot(geneTree, xlab="", sub="",cex=0.3)

minModuleSize = 50

# Module identification using dynamic tree cut
dynamicMods = cutreeDynamic(dendro = geneTree,
                            method="tree", minClusterSize = minModuleSize)

dynamicColors = labels2colors(dynamicMods)
dynamicColors[dynamicColors %in% c("green")] <- "grey"
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", 
                    dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, 
                    guideHang = 0.05, main = "Gene dendrogram and module colors")

### Tune and merge modules
merged.moduleColors <- dynamicColors
merged.moduleColors[merged.moduleColors %in% c("yellow", "turquoise")] <- "yellow"
merged.moduleColors[merged.moduleColors %in% c("green")] <- "grey"

## top module genes
module.genesets <- lapply(all_module_genes, function(r){
  r
}) %>% dplyr::bind_rows()
saveRDS(module.genesets, "./Trajectory/obj/tal.module.genesets.rds")
