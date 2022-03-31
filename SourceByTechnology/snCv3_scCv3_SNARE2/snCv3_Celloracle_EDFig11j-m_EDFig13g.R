# snCv3 - Celloracle (ED Fig 11j-m, ED Fig 13g) ------------------------------------------
library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(patchwork)
library(cicero)
library(SeuratWrappers)
library(monocle3)
set.seed(1234)


####TAL Trajectory Modules
###Prepare seurat objects
load("Kidney_Healthy-Injury_Cell_Atlas_SNARE2-RNA-AC_Seurat_03282022.rda")

# Create a folder to save results
output_folder <- "Celloracle/"

# 2. Load data and make Cell Data Set (CDS) object 
## 2.1. Process data to make CDS object
clusters <- Idents(object = KSAC)
umap.coordinates <- Embeddings(object = KSAC, reduction = "umap.atac")

# prepare matrix data
DefaultAssay(KSAC) <- "ATAC"
count.matrix <- GetAssayData(object = KSAC, slot = "counts")
count.matrix@x[count.matrix@x > 0] <- 1 #binarize counts

#update counts matrix
rownames <- rownames(count.matrix)
chrom <- sapply(rownames, function(x) strsplit(x, split = "-")[[1]][[1]])
loc_start <- sapply(rownames, function(x) strsplit(x, split = "-")[[1]][[2]])
loc_end <- sapply(rownames, function(x) strsplit(x, split = "-")[[1]][[3]])

unique(chrom)

pData <- data.frame(clusters)
fData <- data.frame(site_name = rownames(count.matrix), chromosome = chrom, bp1 = loc_start, bp2 = loc_end)
input_cds <-  suppressWarnings(new_cell_data_set(count.matrix,
                                                 cell_metadata = pData,
                                                 gene_metadata = fData))
input_cds <- monocle3::detect_genes(input_cds)
#Ensure there are no peaks included with zero reads
input_cds <- input_cds[Matrix::rowSums(exprs(input_cds)) != 0,] 


# 3. Quality check and Filtering
#Load conns from "SNARE2 - Combined Seurat object"
conns

# 7. Save results for the next step
all_peaks <- row.names(exprs(input_cds))
all_peaks <- gsub("-","_",all_peaks)
write.csv(x = all_peaks, file = paste0(output_folder, "all_peaks.csv"))
conns$Peak1 <- gsub("-","_",conns$Peak1)
conns$Peak2 <- gsub("-","_",conns$Peak2)
write.csv(x = conns, file = paste0(output_folder, "cicero_connections.csv"))


###Prepare snCv3 RNA Object
KBR <- LoadH5Seurat("Kidney_Healthy-Injury_Cell_Atlas_snCv3_Seurat_03282022.h5Seurat")
Idents(object = KBR) <- "clusters"
Idents(object = KBR) <- factor(Idents(object = KBR), levels = 1:119)
KBR <- subset(KBR, idents = 1:100)

tal.pseudo <- readRDS("tal.psuedotime.rds")
cells.use <- rownames(tal.pseudo)[rownames(tal.pseudo) %in% colnames(KBR)]

KBR.TAL <- subset(KBR, cells = cells.use)
umap.coordinates <- readRDS("tal.umap.embeddings.rds")
umap.coordinates <- umap.coordinates[cells.use,]
dim(umap.coordinates)

KBR.TAL[["umap.traj"]] <- CreateDimReducObject(embeddings = umap.coordinates, key = "umaptraj_", assay = DefaultAssay(KBR.TAL))

tal.mod <- readRDS("tal.sn.sc.module.assignment.rds")
tal.mod <- tal.mod[rownames(KBR.TAL@meta.data)]
KBR.TAL$modules <- tal.mod

tal.pseudo <- tal.pseudo[rownames(KBR.TAL@meta.data),]
KBR.TAL$pseudotime <- tal.pseudo

DimPlot(KBR.TAL, reduction = "umap.traj", pt.size = 0.5, label = TRUE, group.by = "subclass.l3", repel = TRUE) + NoLegend()
DimPlot(KBR.TAL, reduction = "umap.traj", pt.size = 0.5, label = TRUE, group.by = "modules", repel = TRUE) + NoLegend()
FeaturePlot(KBR.TAL, features = "pseudotime", reduction = "umap.traj")

#Find Variable features 
library(pagoda2)
countMatrix <- GetAssayData(object = KBR.TAL, slot = "counts")
countMatrix.batch <- as.factor(unlist(lapply(colnames(countMatrix),function(x) unlist(strsplit(x,"_"))[1])))
names(countMatrix.batch) <- colnames(countMatrix)
p2 <- Pagoda2$new(x = countMatrix, n.cores = 6, trim=10, batch=countMatrix.batch)
p2$adjustVariance(plot = T, gam.k = 10)
#2579 overdispersed genes ... 2579 persisting ... done
var.info <- p2$misc$varinfo; var.info <- var.info[order(var.info$lp, decreasing = F),];
sn.od.genes <- rownames(var.info)[1:1000]

#Load TAL module gene sets (Supplementary Tables 17,20,21)
tal.genes <- read.delim("aTAL_Module_Latent-Time_TFBS_Gene_Set.txt", header = FALSE)
tal.genes <- rownames(KBR.TAL)[rownames(KBR.TAL) %in% tal.genes$V1]
tal.genes <- unique(c(sn.od.genes, tal.genes))

#remove scale.data, subset to var genes
KBR.TAL.sub <- DietSeurat(KBR.TAL, counts = TRUE, 
                          data = TRUE, 
                          scale.data = FALSE,
                          features = tal.genes,
                          dimreducs = "umap.traj")

KBR.TAL.sub@meta.data <- KBR.TAL.sub@meta.data[,c("subclass.l3", "modules", "pseudotime")]

Idents(KBR.TAL.sub) <- "modules"
saveRDS(KBR.TAL.sub, file = paste0(output_folder, "aTAL_Seurat.Rds"))

###Run Celloracle
#See Jupyter Notebook "snCv3_Celoracle_aTAL.ipynb"
###





####STR Subclasses
###Prepare snCv3 RNA Object
Idents(KBR) <- "subclass.l3"
KBR.STR <- subset(KBR, idents = c("FIB", "aFIB", "MYOF"))
KBR.STR <- RunUMAP(object = KBR.STR, reduction = "pca", dims = 1:50, n.neighbors = 30L,
                   min.dist = 0.3)
DimPlot(KBR.STR, reduction = "umap", pt.size = 0.5, label = TRUE, group.by = "subclass.l3", repel = TRUE) + NoLegend()

#Find Variable features 
library(pagoda2)
countMatrix <- GetAssayData(object = KBR.STR, slot = "counts")
countMatrix.batch <- as.factor(unlist(lapply(colnames(countMatrix),function(x) unlist(strsplit(x,"_"))[1])))
names(countMatrix.batch) <- colnames(countMatrix)
p2 <- Pagoda2$new(x = countMatrix, n.cores = 6, trim=10, batch=countMatrix.batch)
p2$adjustVariance(plot = T, gam.k = 10)
#891 overdispersed genes ... 891 persisting ... done
var.info <- p2$misc$varinfo; var.info <- var.info[order(var.info$lp, decreasing = F),];
sn.od.genes <- rownames(var.info)[1:891]

#Load TFs expressed in STR subtypes
load("Kidney_SNARE2_Dual_TF_Markers_Exp_STR_clusters.rda")
tfs <- unique(c(tf.markers.exp$gene, tf.markers.exp$gene2))

#Add in marker genes
subclass.l3.markers.sn #from "snCv3/scCv3 - Marker Genes"
markers <- subclass.l3.markers.sn[subclass.l3.markers.sn$cluster %in% c("FIB","aFIB","MYOF","cycMYOF") & 
                                    subclass.l3.markers.sn$p_val_adj < 0.05,]$gene
all.genes <- unique(c(sn.od.genes,tfs,markers,"SMAD1","SMAD5","SMAD8")) #include SMADS which show low expression

#remove scale.data, subset to var genes
KBR.STR.sub <- DietSeurat(KBR.STR, counts = TRUE, 
                          data = TRUE, 
                          scale.data = FALSE,
                          features = all.genes,
                          dimreducs = "umap")

KBR.STR.sub@meta.data <- KBR.STR.sub@meta.data[,c("subclass.l3","subclass.l2")]
saveRDS(KBR.STR.sub, file = paste0(output_folder, "aSTR_Seurat.Rds"))

###Run Celloracle
#See Jupyter Notebook "snCv3_Celoracle_aSTR.ipynb"
###

