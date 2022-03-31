# snCv3 - Clustering and Annotation -----------------------------------------------------

###Preparation of Individual Count Matrices (example code only)

library(dplyr)
library(Seurat)
library(DropletUtils)
library(Matrix)

#KPMP_20191204A_10X-R	S-1908-000952-R1
dir1 <- "/20191204_KPMP_S-1908-000952-R1_SI-GA-H1/outs/filtered_feature_bc_matrix/"
dir2 <- "/KPMP_10X/RNA_RDS_Files/" 

KB1 <- Read10X(data.dir = dir1)
colnames(KB1) <- paste("KB1", colnames(KB1), sep="_")
saveRDS(KB1, file = paste0(dir2,"KPMP_20191204A_10X-R_UMI_counts.rds"))

#Identify and remove mitochodrial genes 
NK.ms <- KB1
mt.genes <- grep("^MT-", rownames(NK.ms), value = T)
NK.ms.noMT <- NK.ms[!(rownames(NK.ms) %in% mt.genes),]
dim(NK.ms.noMT)
#33525  6759

saveRDS(NK.ms, file = paste0(dir2,"KPMP_20191204A_10X-R_UMI_counts_EmptyBC_Filter.rds"))
saveRDS(NK.ms.noMT, file = paste0(dir2,"KPMP_20191204A_10X-R_UMI_counts_EmptyBC_MT-Filter.rds"))
writeMM(NK.ms.noMT, file= paste0(dir2,"KPMP_20191204A_10X-R_UMI_counts_EmptyBC_MT-Filter.mtx"))

#Run DoubletDetection on Python using the .mtx files

#Remove doublets
countMatrix <- readRDS(paste0(dir2, "KPMP_20191204A_10X-R_UMI_counts_EmptyBC_MT-Filter.rds"))
doublet.det<-read.delim(paste0(dir2, "KPMP_20191204A_10X-R_DoubletDetection_doublets.txt"),sep="\t",header=FALSE)
doublet.det<-as.character(doublet.det$V1)
names(doublet.det)<-colnames(countMatrix)

doublets<-names(doublet.det)[doublet.det == "1"]
singlets<-names(doublet.det)[doublet.det == "0"]

countMatrix.dd<-countMatrix[,!colnames(countMatrix) %in% doublets]
dim(countMatrix.dd)
#33525  6188

saveRDS(countMatrix.dd, file = paste0(dir2,"KPMP_20191204A_10X-R_EmptyBC_MT_DD_Filter.rds"))

#Run for all samples




###Combine all BUKMAP/KPMP filtered counts
KPMPpath = "/KPMP_10X/RNA_RDS_Files"
KPMP_CM.list = list.files(path=KPMPpath, pattern="*_EmptyBC_MT_DD_Filter.rds") 
KPMP_CM.list = paste(KPMPpath,"/",KPMP_CM.list, sep = "")
HuBpath = "BUKMAP_10X/RNA_RDS_Files"
HuB_CM.list = list.files(path=HuBpath, pattern="*_EmptyBC_MT_DD_Filter.rds") 
HuB_CM.list = paste(HuBpath,"/",HuB_CM.list, sep = "")
CM.list = c(KPMP_CM.list,HuB_CM.list)

CM_files_df <- lapply(CM.list, function(x) {readRDS(file = x)})

merge.sparse = function(listMatrixes) {
  # takes a list of sparse matrixes with different columns and adds them row wise
  
  allRownames <- sort(unique(unlist(lapply(listMatrixes,rownames))))
  allColnames <- as.character(unlist(lapply(listMatrixes,colnames)))
  for (currentMatrix in listMatrixes) {
    newRowLocations <- match(rownames(currentMatrix),allRownames)
    indexes <- which(currentMatrix>0, arr.ind = T)
    newRows <- newRowLocations[indexes[,1]]
    columns <- indexes[,2]
    newMatrix <- sparseMatrix(i=newRows,j=columns, x=currentMatrix@x,
                              dims=c(length(allRownames),max(columns)))
    if (!exists("matrixToReturn")) {
      matrixToReturn <- newMatrix
    }
    else {
      matrixToReturn <- cbind2(matrixToReturn,newMatrix)
    }
  }
  rownames(matrixToReturn) <- allRownames
  colnames(matrixToReturn) <- allColnames
  matrixToReturn  
}
KB <- merge.sparse(CM_files_df)




###Create Seurat Object
KBR <- CreateSeuratObject(counts = KB, project = "KPMP Biopsy and Reference", min.cells = 3, min.features = 400)

#Add Metadata
meta <- read.delim("Kidney_KPMP-Biopsy_Ref_10X-R_092020_Meta.txt",sep="\t",header=TRUE,row.names=1)
KBR@meta.data <- meta[rownames(KBR@meta.data),]
head(KBR@meta.data)

#Gene Filter >400 and <7500 genes
KBR <- subset(x = KBR, subset = nFeature_RNA > 400 & nFeature_RNA < 7500)

#Gene Molecule filters
countMatrix <- GetAssayData(object = KBR, slot = "counts")
KB.gmcf <- gene.vs.molecule.cell.filter(countMatrix,min.cell.size=200)
KBR <- subset(KBR, cells = colnames(KB.gmcf))

#Normalize and run PCA
KBR <- NormalizeData(object = KBR, normalization.method = "LogNormalize", 
                     scale.factor = 10000)
KBR <- FindVariableFeatures(object = KBR, selection.method = "vst", nfeatures = 3000)
all.genes <- rownames(x = KBR)
KBR <- ScaleData(object = KBR, features = all.genes)
KBR <- RunPCA(object = KBR, features = VariableFeatures(object = KBR))




###Pagoda2 Clustering
#Filtered Count Matrix from Seurat
countMatrix <- GetAssayData(object = KBR, slot = "counts")

#Batch annotations
countMatrix.batch <- as.factor(unlist(lapply(colnames(countMatrix),function(x) unlist(strsplit(x,"_"))[1])))
names(countMatrix.batch) <- colnames(countMatrix)

#Generate a new pagoda2 object
p2 <- Pagoda2$new(x = countMatrix, n.cores = 6, trim=10, batch=countMatrix.batch)

# Adjust the variance
p2$adjustVariance(plot = T, gam.k = 10)
#5526 overdispersed genes ... 5526 persisting ... done

#Calculate a PCA reduction with the number of PCs specified by nPCs
# and using only the n.odgenes overdispersed genes 
p2$calculatePcaReduction(nPcs = 50, n.odgenes = 5526, maxit=1000)

#Generate K-nearest neighbour graph
p2$makeKnnGraph(k = 50, type = 'PCA', center = T, weight.type = 'none', n.cores = 6, distance = 'cosine')

# Identify clusters using the infomap.community method
# on the basis of the reduction called 'PCA' (generated above)
# Save the resulting clustering as 'infomap'
p2$getKnnClusters(method = igraph::infomap.community, type = 'PCA', name = 'infomap')

saveRDS(p2,file = "Kidney_KPMP-Biopsy_Ref_10X-R_092020_Pagoda_50PC_all-var-genes_k50.RDS")
#saveRDS(p2,file = "Kidney_KPMP-Biopsy_Ref_10X-R_092020_Pagoda_50PC_k100_all-var-genes.RDS")
#saveRDS(p2,file = "Kidney_KPMP-Biopsy_Ref_10X-R_092020_Pagoda_50PC_k200_all-var-genes.RDS")
#saveRDS(p2,file = "Kidney_KPMP-Biopsy_Ref_10X-R_092020_Pagoda_50PC_k500_all-var-genes.RDS")

k50infomap <- p2$clusters$PCA$infomap
#k100infomap <- p2$clusters$PCA$infomap
#k200infomap <- p2$clusters$PCA$infomap
#k500infomap <- p2$clusters$PCA$infomap

#pagoda2 clusters and PCA values for umap
p2 <- readRDS("Kidney_KPMP-Biopsy_Ref_10X-R_092020_Pagoda_50PC_k100_all-var-genes.RDS")
cell.embeddings <- p2$reductions$PCA[rownames(Embeddings(KBR, reduction = "pca")),]
KBR[["pca"]] <- CreateDimReducObject(embeddings = cell.embeddings, key = "pca_", assay = DefaultAssay(KBR))

KBR <- RunUMAP(object = KBR, reduction = "pca", dims = 1:50, n.neighbors = 30L,
               min.dist = 0.3)

KBR[["pagoda_k50_infomap"]] <- k50infomap[rownames(KBR@meta.data)]
KBR[["pagoda_k100_infomap"]] <- k100infomap[rownames(KBR@meta.data)]
KBR[["pagoda_k200_infomap"]] <- k200infomap[rownames(KBR@meta.data)]
KBR[["pagoda_k500_infomap"]] <- k500infomap[rownames(KBR@meta.data)]



###Find marker gene expression at different k values (Identify top 5 distinct markers)
Idents(KBR) <- "pagoda_k100_infomap"
Idents(KBR) <- factor(Idents(KBR), levels = 1:122)
KBR.markers <- FindAllMarkers(KBR, idents = 1:122, only.pos = TRUE, max.cells.per.ident = 1000,
                              logfc.threshold = 0.25, min.pct = 0.25)
cl.mark <- KBR.markers[KBR.markers$p_val_adj < 0.05,]
cl.mark <- cl.mark[cl.mark$pct.2 < 0.1,]
cl.mark <- cl.mark[cl.mark$pct.1 > 0.3,]
cl.mark <- cl.mark[cl.mark$avg_logFC > 1,]
cl.mark <- distinct(cl.mark, gene, .keep_all = TRUE)
cl.mark %>% group_by(cluster) %>% top_n(5, avg_logFC) -> top5




###Add in Ribosomal gene proportions
ER.GO <- read.delim("Cytosolic_Ribosome_Gene_Ontology.txt",sep="\t",header=TRUE,row.names=1)
ER.GO <- as.character(ER.GO$Symbol)
ER.GO <- intersect(ER.GO, rownames(KBR))
KBR[["percent.er"]] <- PercentageFeatureSet(KBR, features = ER.GO)
VlnPlot(KBR, features = c("nFeature_RNA", "nCount_RNA", "percent.er"), 
        ncol = 1, pt.size = -1)

#Add in mitochondrial transcript percentages (using pre-MT filtered count matrices)
dir <- "/KPMP_10X/RNA_RDS_Files/" 
KPMPpath = "/KPMP_10X/RNA_RDS_Files" 
KPMP_CM.list = list.files(path=KPMPpath, pattern="*UMI_counts_EmptyBC_Filter.rds") 
KPMP_CM.list = paste(KPMPpath,"/",KPMP_CM.list, sep = "")
HuBpath = "/BUKMAP_10X/RNA_RDS_Files"
HuB_CM.list = list.files(path=HuBpath, pattern="*UMI_counts_EmptyBC_Filter.rds") 
HuB_CM.list = paste(HuBpath,"/",HuB_CM.list, sep = "")
CM.list = c(KPMP_CM.list,HuB_CM.list)
CM_files_df <- lapply(CM.list, function(x) {readRDS(file = x)})
merge.sparse = function(listMatrixes) {
  # takes a list of sparse matrixes with different columns and adds them row wise
  
  allRownames <- sort(unique(unlist(lapply(listMatrixes,rownames))))
  allColnames <- as.character(unlist(lapply(listMatrixes,colnames)))
  for (currentMatrix in listMatrixes) {
    newRowLocations <- match(rownames(currentMatrix),allRownames)
    indexes <- which(currentMatrix>0, arr.ind = T)
    newRows <- newRowLocations[indexes[,1]]
    columns <- indexes[,2]
    newMatrix <- sparseMatrix(i=newRows,j=columns, x=currentMatrix@x,
                              dims=c(length(allRownames),max(columns)))
    if (!exists("matrixToReturn")) {
      matrixToReturn <- newMatrix
    }
    else {
      matrixToReturn <- cbind2(matrixToReturn,newMatrix)
    }
  }
  rownames(matrixToReturn) <- allRownames
  colnames(matrixToReturn) <- allColnames
  matrixToReturn  
}
KB <- merge.sparse(CM_files_df)
KB <- CreateSeuratObject(counts = KB, project = "Kidney Atlas", min.cells = 3, min.features = 400)
KB[["percent.mt"]] <- PercentageFeatureSet(KB, pattern = "^MT-")
meta <- KB@meta.data
KBR$percent.mt <- meta[rownames(KBR@meta.data),]$percent.mt

VlnPlot(KBR, features = c("nFeature_RNA", "percent.er", "percent.mt"), 
        ncol = 1, pt.size = -1)


###Annotations and low qualtity trimming
##See Methods section for full description on cell type annotations 

###Identify REN
ren <- WhichCells(object = KBR, expression = REN > 3)
Idents(KBR) <- "class"
stroma <- WhichCells(object = KBR, idents = "stroma cells")
ren <- ren[ren %in% stroma]


###Separation of AEA from DVR (subclustering)
library(pagoda2)
require(parallel)
KBR.aea <- subset(KBR, idents = 63) #EC-AEA/DVR annotated cluster - different from final cluster #
countMatrix <- GetAssayData(object = KBR.aea, slot = "counts")
countMatrix.batch <- as.factor(unlist(lapply(colnames(countMatrix),function(x) unlist(strsplit(x,"_"))[1])))
names(countMatrix.batch) <- colnames(countMatrix)
p2 <- Pagoda2$new(x = countMatrix, n.cores = 6, trim=10, batch=countMatrix.batch)
p2$adjustVariance(plot = T, gam.k = 10)
#503 overdispersed genes ... 503 persisting ... done
p2$calculatePcaReduction(nPcs = 10, n.odgenes = 500, maxit=1000)
p2$makeKnnGraph(k = 100, type = 'PCA', center = T, weight.type = 'none', n.cores = 6, distance = 'cosine')
p2$getKnnClusters(method = igraph::infomap.community, type = 'PCA', name = 'infomap')
k100infomap <- p2$clusters$PCA$infomap
KBR.aea[["pagoda_k100_infomap_AEA"]] <- k100infomap[rownames(KBR.aea@meta.data)]
DimPlot(KBR.aea, reduction = "umap", pt.size = 0.5, label = TRUE, group.by = "pagoda_k100_infomap_AEA", repel = TRUE) + NoLegend()
Idents(KBR.aea) <- "pagoda_k100_infomap_AEA"
CL.AEA <- WhichCells(KBR.aea, idents = c(5,7))
CL.DVR <- WhichCells(KBR.aea, idents = c(5,7), invert = TRUE)


###
#Final Seurat object (with updated metadata): Kidney_Healthy-Injury_Cell_Atlas_snCv3_Seurat_03282022.h5Seurat
###

# snCv3/SNARE2 - Data Integration -----------------------------------------
library(Seurat)
library(SeuratDisk)
library(patchwork)
library(dplyr)
library(Matrix)
library(ggplot2)

###Prepare Ref Data
KBR <- LoadH5Seurat("Kidney_Healthy-Injury_Cell_Atlas_snCv3_Seurat_03282022.h5Seurat")

p2 #pagoda2 object generated when clustering snCv3 data (not provided)
var.info <- p2$misc$varinfo; var.info <- var.info[order(var.info$lp, decreasing = F),];
sn.od.genes <- rownames(var.info)[1:5526]

KBR <- SCTransform(KBR, verbose = TRUE)
VariableFeatures(KBR) <- sn.od.genes
KBR <- RunPCA(KBR, verbose = FALSE)
pca.loadings <- p2$misc$PCA$v[sn.od.genes,]
colnames(pca.loadings) <- paste("PC_", 1:50, sep = "")
KBR@reductions$pca@feature.loadings <- pca.loadings
cell.embeddings <- p2$reductions$PCA[rownames(Embeddings(KBR, reduction = "pca")),]
colnames(cell.embeddings) <- paste("PC_", 1:50, sep = "")
KBR[["pca"]] <- CreateDimReducObject(embeddings = cell.embeddings, key = "pca_", assay = DefaultAssay(KBR))
KBR <- RunUMAP(object = KBR, reduction = "pca", dims = 1:50, n.neighbors = 30L,
               min.dist = 0.3, umap.method = "uwot-learn")
KBR <- FindNeighbors(
  object = KBR,
  reduction = "pca",
  dims = 1:50,
  graph.name = "pca.annoy.neighbors", 
  k.param = 100,
  cache.index = TRUE,
  return.neighbor = TRUE,
  l2.norm = TRUE
)

#FindNeighbors wipes the pca loadings, so re-add it
KBR@reductions$pca@feature.loadings <- pca.loadings
Loadings(object = KBR[["pca"]])




###Prepare Query Data
##Load in all SNARE RNA count matrices
dir1 = "/BUKMAP_SNARE/RNA_RDS_Files/"
KID.list = list.files(path=dir1, pattern = "_dual.rds$", include.dirs= TRUE)
CM_files_df <- lapply(KID.list, function(x) {readRDS(file = x)})

#Combine all Samples
merge.sparse = function(listMatrixes) {
  # takes a list of sparse matrixes with different columns and adds them row wise
  
  allRownames <- sort(unique(unlist(lapply(listMatrixes,rownames))))
  allColnames <- as.character(unlist(lapply(listMatrixes,colnames)))
  for (currentMatrix in listMatrixes) {
    newRowLocations <- match(rownames(currentMatrix),allRownames)
    indexes <- which(currentMatrix>0, arr.ind = T)
    newRows <- newRowLocations[indexes[,1]]
    columns <- indexes[,2]
    newMatrix <- sparseMatrix(i=newRows,j=columns, x=currentMatrix@x,
                              dims=c(length(allRownames),max(columns)))
    if (!exists("matrixToReturn")) {
      matrixToReturn <- newMatrix
    }
    else {
      matrixToReturn <- cbind2(matrixToReturn,newMatrix)
    }
  }
  rownames(matrixToReturn) <- allRownames
  colnames(matrixToReturn) <- allColnames
  matrixToReturn  
}
KB <- merge.sparse(CM_files_df)

#Gene Molecule filter
KID.SNARE <- CreateSeuratObject(counts = KB, project = "Kidney SNARE2", min.cells = 3, min.features = 200)
countMatrix <- GetAssayData(object = KID.SNARE, slot = "counts")
KID.SNARE.gmcf <- gene.vs.molecule.cell.filter(countMatrix,min.cell.size=200)
KID.SNARE <- subset(KID.SNARE, cells = colnames(KID.SNARE.gmcf))
KID.SNARE <- SCTransform(KID.SNARE, verbose = TRUE)



###Find anchor genes
anchors <- FindTransferAnchors(
  reference = KBR,
  query = KID.SNARE,
  normalization.method = "SCT",
  k.filter = NA,
  reference.reduction = "pca",
  reference.neighbors = "pca.annoy.neighbors",
  dims = 1:50
)

KID.SNARE <- MapQuery(
  anchorset = anchors,
  query = KID.SNARE,
  reference = KBR,
  refdata = list(
    class = "class",
    subclass.l1 = "subclass.l1",
    subclass.l2 = "subclass.l2",
    clusters = "clusters"
  ),
  reference.reduction = "pca", 
  reduction.model = "umap"
)


###merge reference and query
KBR$id <- 'reference'
KID.SNARE$id <- 'query'
refquery <- merge(KBR, KID.SNARE)
refquery[["pca"]] <- merge(KBR[["pca"]], KID.SNARE[["ref.pca"]])
refquery <- RunUMAP(object = refquery, reduction = "pca", dims = 1:50, n.neighbors = 30L,
                    min.dist = 0.2)




###Clustering in P2 using integrated PCs
#Filtered Count Matrix from Seurat
countMatrix <- GetAssayData(object = refquery, slot = "counts")

#Batch annotations
countMatrix.batch <- as.factor(unlist(lapply(colnames(countMatrix),function(x) unlist(strsplit(x,"_"))[1])))
names(countMatrix.batch) <- colnames(countMatrix)

#Cluster
p2 <- Pagoda2$new(x = countMatrix, n.cores = 6, trim=10, batch=countMatrix.batch)
p2$adjustVariance(plot = T, gam.k = 10)
#4846 overdispersed genes ... 4846 persisting ... done
p2$calculatePcaReduction(nPcs = 50, n.odgenes = 4846, maxit=1000)
p2$reductions$PCA <- Embeddings(refquery, reduction = "pca")
p2$makeKnnGraph(k = 100, type = 'PCA', center = T, weight.type = 'none', n.cores = 6, distance = 'cosine')
p2$getKnnClusters(method = igraph::infomap.community, type = 'PCA', name = 'infomap')

#Add pagoda2 clusters 
k100infomap <- p2$clusters$PCA$infomap
refquery[["pagoda_k100_infomap_coembed"]] <- k100infomap[rownames(refquery@meta.data)]

#Add in k = 200 resolution clusters
p2$makeKnnGraph(k = 200, type = 'PCA', center = T, weight.type = 'none', n.cores = 6, distance = 'cosine')
p2$getKnnClusters(method = igraph::infomap.community, type = 'PCA', name = 'infomap')
k200infomap <- p2$clusters$PCA$infomap
refquery[["pagoda_k200_infomap_coembed"]] <- k200infomap[rownames(refquery@meta.data)]






#The SNARE2 (RNA) component of the integrated
#clusters were then annotated to the most overlapping or correlated snCv3 cluster, with manual
#inspection of cell type markers used to confirm identities.  

#Updated metadata tables to include snCv3/scCv3 subclass.l3,subclass.l2,subclass.l1,class,substructure annotations



###Identify REN for query (SNARE)
ren <- WhichCells(object = refquery, expression = REN > 3)
Idents(refquery) <- "class"
stroma <- WhichCells(object = refquery, idents = "stroma cells")
ren <- ren[ren %in% stroma]
Idents(refquery) <- "id"
query <- WhichCells(object = refquery, idents = "query")
ren <- ren[ren %in% query]


###Separate AEA from DVR
library(pagoda2)
require(parallel)
Idents(refquery) <- "clusters"
KBR.aea <- subset(refquery, idents = c(63:64))
countMatrix <- GetAssayData(object = KBR.aea, slot = "counts")
countMatrix.batch <- as.factor(unlist(lapply(colnames(countMatrix),function(x) unlist(strsplit(x,"_"))[1])))
names(countMatrix.batch) <- colnames(countMatrix)
p2 <- Pagoda2$new(x = countMatrix, n.cores = 6, trim=10, batch=countMatrix.batch)
p2$adjustVariance(plot = T, gam.k = 10)
#661 overdispersed genes ... 661 persisting ... done
p2$calculatePcaReduction(nPcs = 10, n.odgenes = 661, maxit=1000)
p2$makeKnnGraph(k = 100, type = 'PCA', center = T, weight.type = 'none', n.cores = 6, distance = 'cosine')
p2$getKnnClusters(method = igraph::infomap.community, type = 'PCA', name = 'infomap')
k100infomap <- p2$clusters$PCA$infomap
KBR.aea[["pagoda_k100_infomap_AEA"]] <- k100infomap[rownames(KBR.aea@meta.data)]
KBR.aea <- RunUMAP(object = KBR.aea, reduction = "pca", dims = 1:50, n.neighbors = 30L,
                   min.dist = 0.3)
DimPlot(KBR.aea, reduction = "umap", pt.size = 0.5, label = TRUE, group.by = "pagoda_k100_infomap_AEA", repel = TRUE) + NoLegend()
DimPlot(KBR.aea, reduction = "umap", pt.size = 0.5, label = TRUE, group.by = "region.l2", repel = TRUE) + NoLegend()
DimPlot(KBR.aea, reduction = "umap", pt.size = 0.5, label = TRUE, group.by = "subclass.l3", repel = TRUE) + NoLegend()
Idents(KBR.aea) <- "pagoda_k100_infomap_AEA"
#Annotate pagoda_k100_infomap_AEA clusters based on overlapping snCv3 annotations to either AEA or DVR





###Remove low quality, ambiguous data and doublets
dir1 = "Doublets/" #From snarePip

doublet.list = list.files(path=dir1, pattern = ".rds$", include.dirs= TRUE)
setwd(dir1)

meta <- KID.SNARE@meta.data

getDoubletMeta <- function(meta, doublet.list) {
  meta <- do.call(rbind, lapply(1:length(doublet.list), function(ct) {
    print(paste("Running for exp:", ct))
    doublet <- readRDS(doublet.list[[ct]])
    cells <- doublet$doubletDetection$cells
    cells.keep <- cells[cells %in% rownames(meta)]
    meta.sub <- meta[cells,]
    meta.sub$doublet.detection <- doublet$doubletDetection$V2
    cells <- names(doublet$scrublet$is_doublet)
    meta.sub <- meta.sub[cells,]
    meta.sub$scrublet <- doublet$scrublet$is_doublet
    meta.sub <- meta.sub[cells.keep,]
    meta.sub
  }))
  return(meta)
}

meta.new <- getDoubletMeta(meta, doublet.list)
meta.new <- meta.new[rownames(meta),]

refquery@meta.data$scrublet <- rep("FALSE", length(rownames(refquery@meta.data)))
refquery@meta.data$doublet.detection <- rep(0, length(rownames(refquery@meta.data)))

ref.cells <- rownames(refquery@meta.data[refquery@meta.data$id == "query",])
refquery@meta.data[ref.cells,]$scrublet <- meta.new[ref.cells,]$scrublet
refquery@meta.data[ref.cells,]$doublet.detection <- meta.new[ref.cells,]$doublet.detection

#Remove ambiguous or low quality cells
refquery <- subset(refquery, idents = 1:100)

#Remove doublets identified by DoubletDetection and scrublet (more stringent doublet identification)
doublets <- rownames(refquery@meta.data[refquery@meta.data$doublet.detection == 1 & 
                                          refquery@meta.data$scrublet == "TRUE",])
refquery <- subset(refquery, cells = doublets, invert = TRUE)



###


# snCv3/scCv3 - Data Integration ------------------------------------------
library(Seurat)
library(SeuratDisk)
library(patchwork)
library(dplyr)
library(Matrix)
library(ggplot2)
library(pagoda2)
library(igraph)


###Prepare Ref Data
KBR <- LoadH5Seurat("Kidney_Healthy-Injury_Cell_Atlas_snCv3_Seurat_03282022.h5Seurat")

p2 #pagoda2 object generated when clustering snCv3 data (not provided)
var.info <- p2$misc$varinfo; var.info <- var.info[order(var.info$lp, decreasing = F),];
sn.od.genes <- rownames(var.info)[1:5526]


KBR <- SCTransform(KBR, verbose = TRUE)
VariableFeatures(KBR) <- sn.od.genes
KBR <- RunPCA(KBR, verbose = FALSE)
pca.loadings <- p2$misc$PCA$v[sn.od.genes,]
colnames(pca.loadings) <- paste("PC_", 1:50, sep = "")
KBR@reductions$pca@feature.loadings <- pca.loadings
cell.embeddings <- p2$reductions$PCA[rownames(Embeddings(KBR, reduction = "pca")),]
colnames(cell.embeddings) <- paste("PC_", 1:50, sep = "")
KBR[["pca"]] <- CreateDimReducObject(embeddings = cell.embeddings, key = "pca_", assay = DefaultAssay(KBR))
KBR <- RunUMAP(object = KBR, reduction = "pca", dims = 1:50, n.neighbors = 30L,
               min.dist = 0.3, umap.method = "uwot-learn")
KBR <- FindNeighbors(
  object = KBR,
  reduction = "pca",
  dims = 1:50,
  graph.name = "pca.annoy.neighbors", 
  k.param = 100,
  cache.index = TRUE,
  return.neighbor = TRUE,
  l2.norm = TRUE
)

#FindNeighbors wipes the pca loadings, so re-add it
KBR@reductions$pca@feature.loadings <- pca.loadings
Loadings(object = KBR[["pca"]])

KBR.SCT <- KBR
rm(KBR)



###Prepare Query Data
sc.kpmp <- readRDS("Premiere_LD_RawCounts.RDS")  #Raw counts for scCv3, not provided
sckb <- CreateSeuratObject(counts = sc.kpmp, project = "sc KPMP Biopsy", min.cells = 3, min.features = 400)#,
sckb <- SCTransform(sckb, verbose = FALSE)

#use pagoda2 to find variable features
countMatrix <- GetAssayData(object = sckb, slot = "counts")
countMatrix.batch <- as.factor(unlist(lapply(colnames(countMatrix),function(x) unlist(strsplit(x,"_"))[1])))
names(countMatrix.batch) <- colnames(countMatrix)
p2.sc <- Pagoda2$new(x = countMatrix, n.cores = 6, trim=10, batch = countMatrix.batch)
p2.sc$adjustVariance(plot = T, gam.k = 10)
#3977 overdispersed genes ... 3977 persisting ... done
var.info <- p2.sc$misc$varinfo; var.info <- var.info[order(var.info$lp, decreasing = F),];
sc.od.genes <- rownames(var.info)[1:3977]
VariableFeatures(sckb) <- sc.od.genes

#Run PCA and UMAP
sckb <- RunPCA(sckb, features = VariableFeatures(object = sckb))
sckb <- RunUMAP(sckb, dims = 1:50)


###Find anchor genes
anchors <- FindTransferAnchors(
  reference = KBR.SCT,
  query = sckb,
  normalization.method = "SCT",
  k.filter = NA,
  reference.reduction = "pca",
  reference.neighbors = "pca.annoy.neighbors",
  dims = 1:50
)

sckb <- MapQuery(
  anchorset = anchors,
  query = sckb,
  reference = KBR.SCT,
  refdata = list(
    class = "class",
    subclass.l1 = "subclass.l1",
    subclass.l2 = "subclass.l2",
    subclass.l3 = "subclass.l3"
  ),
  reference.reduction = "pca", 
  reduction.model = "umap"
)


###merge reference and query
KBR.SCT$id <- 'reference'
sckb$id <- 'query'
refquery <- merge(KBR.SCT, sckb)
refquery[["pca"]] <- merge(KBR.SCT[["pca"]], sckb[["ref.pca"]])
refquery <- RunUMAP(object = refquery, reduction = "pca", dims = 1:50, n.neighbors = 30L,
                    min.dist = 0.2)
umap.coordinates <- rbind(Embeddings(object = KBR.SCT, reduction = "umap"),
                          Embeddings(object = sckb, reduction = "ref.umap"))
refquery[["ref.umap"]] <- CreateDimReducObject(embeddings = umap.coordinates, key = "refumap_", assay = DefaultAssay(refquery))


###Subset out ambiguous or low quality clusters from snCv3 data
to.remove <- rownames(refquery@meta.data[refquery@meta.data$clusters %in% c(101:119),])
refquery <- subset(refquery, cells = to.remove, invert = TRUE)


###Clustering in P2 using integrated PCs
#Filtered Count Matrix from Seurat
countMatrix <- GetAssayData(object = refquery, slot = "counts")

#Generate a new pagoda2 object
p2 <- Pagoda2$new(x = countMatrix, n.cores = 6, trim=10)

#Adjust the variance
p2$adjustVariance(plot = T, gam.k = 10)
#4640 overdispersed genes ... 4640 persisting ... done

#Add in integrated PCA embeddings
p2$reductions$PCA <- Embeddings(refquery, reduction = "pca")

#Generate K-nearest neighbour graph
p2$makeKnnGraph(k = 100, type = 'PCA', center = T, weight.type = 'none', n.cores = 6, distance = 'cosine')

#Identify clusters using the infomap.community method
#on the basis of the reduction called 'PCA' (generated above)
#Save the resulting clustering as 'infomap'
p2$getKnnClusters(method = igraph::infomap.community, type = 'PCA', name = 'infomap')

#Add pagoda2 clusters 
k100infomap <- p2$clusters$PCA$infomap
refquery[["pagoda_k100_infomap_coembed"]] <- k100infomap[rownames(refquery@meta.data)]


#The scCv3 component of the integrated clusters was then annotated to
# the most overlapping or correlated snCv3 subclass, with manual
# inspection of cell type markers used to confirm identities. Cell types that could 
# not be accurately resolved (PT-S1/PT-S2) were kept merged. Integrated clusters that
# overlapped different classes of cell types or that were too ambiguous to annotate were
# considered low quality and were removed from the analysis. 

#Updated metadata tables to include snCv3/scCv3 subclass.l3,subclass.l2,subclass.l1,class,substructure annotations


###Identify REN for scCv3
ren <- WhichCells(object = refquery, expression = REN > 3)
Idents(refquery) <- "class"
stroma <- WhichCells(object = refquery, idents = "stroma cells")
ren <- ren[ren %in% stroma]
Idents(refquery) <- "subclass.l2"
vsmc <- WhichCells(object = refquery, idents = "VSMC")
ren <- ren[ren %in% vsmc]

Idents(refquery) <- "id"
query <- WhichCells(object = refquery, idents = "query")
ren <- ren[ren %in% query]


###Separate AEA from DVR
library(pagoda2)
require(parallel)
Idents(refquery) <- "subclass.l3"
KBR.aea <- subset(refquery, idents = c("EC-AEA/DVR"))
query.KBR.aea <- subset(KBR.aea, cells = rownames(KBR.aea@meta.data[KBR.aea@meta.data$id == "query",]))
countMatrix <- GetAssayData(object = query.KBR.aea, slot = "counts")
countMatrix.batch <- as.factor(unlist(lapply(colnames(countMatrix),function(x) unlist(strsplit(x,"_"))[1])))
names(countMatrix.batch) <- colnames(countMatrix)
p2 <- Pagoda2$new(x = countMatrix, n.cores = 6, trim=10, batch=countMatrix.batch)
p2$adjustVariance(plot = T, gam.k = 10)
#100 overdispersed genes ... 100 persisting ... done
p2$calculatePcaReduction(nPcs = 10, n.odgenes = 100, maxit=1000)
p2$makeKnnGraph(k = 300, type = 'PCA', center = T, weight.type = 'none', n.cores = 6, distance = 'cosine')
p2$getKnnClusters(method = igraph::infomap.community, type = 'PCA', name = 'infomap')
k300infomap <- p2$clusters$PCA$infomap
cell.embeddings <- p2$reductions$PCA
query.KBR.aea[["pca"]] <- CreateDimReducObject(embeddings = cell.embeddings, key = "pca_", assay = DefaultAssay(query.KBR.aea))
query.KBR.aea[["pagoda_k300_infomap_AEA"]] <- k300infomap[rownames(query.KBR.aea@meta.data)]
query.KBR.aea <- RunUMAP(object = query.KBR.aea, reduction = "pca", dims = 1:10, n.neighbors = 30L,
                         min.dist = 0.3)
Idents(query.KBR.aea) <- "pagoda_k300_infomap_AEA"
#Annotate pagoda_k300_infomap_AEA clusters based on correlation with snCv3 annotations to either AEA or DVR



# SNARE2 - Combined Seurat object ----------------------------------------------------------------

library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(dplyr)
library(Matrix)
library(ggplot2)
library(EnsDb.Hsapiens.v86)
set.seed(1234)

fragments <- "BUKMAP.fragments.sort.tsv.gz"

###Prepare SNARE2 RNA
refquery <- LoadH5Seurat("Kidney_Healthy-Injury_Cell_Atlas_Integrated-snCv3-SNARE-RNA_Seurat_08032021.h5Seurat")
snare.cells <- rownames(refquery@meta.data[refquery@meta.data$id == "query",])
KID.SNARE <- subset(refquery, cells = snare.cells)
rna <- GetAssayData(KID.SNARE, assay = "RNA", slot = "counts")
KID.SNARE.AC <- CreateSeuratObject(counts = rna)
meta <- KID.SNARE@meta.data

###AC counts
counts <- readRDS("Kidney_SNARE2_Dual_11302020c.pmat.rds") #Combined peaks across samples from snarePip
counts <- counts@assays@data[[1]]
colnames(counts) <- sub(".*#", "", colnames(counts))

#create and pre-process atac seurat 
KID.SNARE.AC <- subset(KID.SNARE.AC, cells = colnames(counts))
KID.SNARE.AC[['ATAC']] <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = "hg38",
  fragments = fragments
)

###Add Metadata from SNARE-RNA 
KID.SNARE.AC@meta.data <- meta[rownames(KID.SNARE.AC@meta.data),]


###Call Peaks based on RNA cell types and re-generate seurat object
DefaultAssay(KID.SNARE.AC) <- "ATAC"
peaks <- CallPeaks(
  object = KID.SNARE.AC,
  group.by = "clusters",
  macs2.path = "/anaconda3/bin/macs3"
)
peaks.l1 <- CallPeaks(
  object = KID.SNARE.AC,
  group.by = "subclass.l1",
  macs2.path = "/anaconda3/bin/macs3"
)
peaks.l3 <- CallPeaks(
  object = KID.SNARE.AC,
  group.by = "subclass.l3",
  macs2.path = "/anaconda3/bin/macs3"
)

# Create a unified set of peaks to quantify in each dataset
combined.peaks <- reduce(x = c(peaks, peaks.l1, peaks.l3))

# Filter out bad peaks based on length
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
combined.peaks

counts <- FeatureMatrix(fragments = Fragments(KID.SNARE.AC),
                        cells = colnames(KID.SNARE.AC),
                        features = peaks,
                        process_n = 20000,
                        sep = c('-','-'))

# keep only peaks from chr 1-22 or chrX or chrY
chrom <- sapply(rownames(counts), function(x) strsplit(x, split = "-")[[1]][[1]])
unique(chrom)
chrom <- sapply(chrom, function(x) strsplit(x, split = "_")[[1]][[1]])
loc1 <- sapply(rownames(counts), function(x) strsplit(x, split = "-")[[1]][[2]])
loc2 <- sapply(rownames(counts), function(x) strsplit(x, split = "-")[[1]][[3]])
rownames(counts) <- paste(chrom, loc1, loc2, sep = "-")
chrom.selected <- as.factor(chrom)
chrom.remove <- names(chrom.selected[chrom.selected %in% c("chrUn", "chrM")])
chrom.keep <- setdiff(names(chrom.selected), chrom.remove)
filtered.counts<- counts[which(rownames(counts) %in% chrom.keep),]

counts <- filtered.counts


###Re-create RNA/AC seurat object

###RNA counts
KSAC <- CreateSeuratObject(counts = rna, meta.data = meta, assay = "RNA")

DefaultAssay(KSAC) <- "RNA"
KSAC <- NormalizeData(KSAC)
all.genes <- rownames(x = KSAC)
KSAC <- ScaleData(object = KSAC, features = all.genes)

umap.coordinates <- Embeddings(object = KID.SNARE.AC, reduction = "umap.rna")[colnames(KSAC),]
KSAC[["umap.rna"]] <- CreateDimReducObject(embeddings = umap.coordinates, key = "umap.rna_", assay = DefaultAssay(KSAC))


###AC counts
#create and pre-process atac seurat 
KSAC[['ATAC']] <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = "hg38",
  fragments = fragments
)

# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)

# change to UCSC style since the data was mapped to hg19
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

# add the gene information to the object
Annotation(KSAC[["ATAC"]]) <- annotations



###Quality Control
DefaultAssay(KSAC) <- "ATAC"
KSAC <- TSSEnrichment(KSAC)
KSAC <- NucleosomeSignal(KSAC)
KSAC$blacklist_fraction <- FractionCountsInRegion(
  object = KSAC,
  assay = 'ATAC',
  regions = blacklist_hg38
)


#preprocessing
DefaultAssay(KSAC) <- "ATAC"
KSAC <- FindTopFeatures(KSAC, min.cutoff = 10)
KSAC <- RunTFIDF(KSAC) #latent semantic indexing
KSAC <- RunSVD(KSAC)


###Add in motif matrix
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)

# Get a list of motif position frequency matrices from the JASPAR database
file = "JASPAR2020_CORE_vertebrates_non-redundant_pfms_jaspar.txt"
pfm <- readJASPARMatrix(file, matrixClass="PFM")

# add motif information
KSAC <- AddMotifs(
  object = KSAC,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm
)

###Add in Chromvar TF activites
KSAC <- RunChromVAR(
  object = KSAC,
  genome = BSgenome.Hsapiens.UCSC.hg38
)



###Cluster ATAC using cisTopic
library(cisTopic)
DefaultAssay(KSAC) <- "ATAC"
counts <- GetAssayData(object = KSAC, slot = "counts")
rownames(counts) <- sub("-",":",rownames(counts))


# initialize cisTopic object from count matrix
cisTopicObject <- createcisTopicObject(counts, min.cells = 50, min.regions = 200, keepCountsMatrix = FALSE)

# run LDA model
cisTopicObject <- runCGSModels(cisTopicObject, topic=c(2, 10, 20, 30, 40, 50), seed=987, nCores=6, 
                               burnin = 250, iterations = 500, addModels=FALSE)

# select for model
cisTopicObject <- cisTopic::selectModel(cisTopicObject, type = "maximum")
#50

# check likelihood stablization
logLikelihoodByIter(cisTopicObject, select=c(2, 10, 20, 30, 40, 50))

# interpretation
# identification of cell state
cisTopicObject <- runUmap(cisTopicObject, target="cell")

# pull out umap coordinates
umap.coordinates <- cisTopicObject@dr$cell[["Umap"]]

# (Explicitly running umap to get the significance matrix)
modelMat <- modelMatSelection(cisTopicObject, target="cell", method = "Z-score")
set.seed(123)
Umap <- umap::umap(t(modelMat))

# Add Cistopic umap coordinates to Seurat
KSAC <- subset(KSAC, cells = rownames(umap.coordinates))
KSAC[["umap.atac"]] <- CreateDimReducObject(embeddings = umap.coordinates, key = "atacUMAP_", assay = DefaultAssay(KSAC))
KSAC[["topic"]] <- CreateDimReducObject(embeddings = Umap$data, key = "topic_", assay = DefaultAssay(KSAC))





###Create gene activity matrix
gene.activities <- GeneActivity(KSAC)

# add the gene activity matrix to the Seurat object as a new assay and normalize it
KSAC[['Activity']] <- CreateAssayObject(counts = gene.activities)
KSAC <- NormalizeData(
  object = KSAC,
  assay = 'Activity',
  normalization.method = 'LogNormalize',
  scale.factor = median(KSAC$nCount_RNA)
)
DefaultAssay(KSAC) <- 'Activity'
KSAC <- ScaleData(KSAC)

Idents(KSAC) <- "clusters"
Idents(KSAC) <- factor(Idents(KSAC), levels = 1:100)




###Add in cicero connections
library(cicero)
library(SeuratWrappers)
library(monocle3)

KID.cds <- as.cell_data_set(x = KSAC)
KID.cicero <- make_cicero_cds(KID.cds, reduced_coordinates = reducedDims(KID.cds)$UMAP.ATAC)

#Cicero analysis: generate predicted gene activity matrix from chrom data
DefaultAssay(KSAC) <- 'ATAC'
Idents(KSAC) <- "clusters"
Idents(KSAC) <- factor(Idents(object = KSAC), levels = 1:100)

## Run cicero
hg38.chr.lengths <- read.table("hg38.chr.lengths.txt", header = F, sep = "\t") ## A tab separated text file with chromosome lengths
hg38.chr.lengths[[2]] <- as.numeric(hg38.chr.lengths[[2]])

clusters <- Idents(object = KSAC)
count.matrix <- GetAssayData(object = KSAC, slot = "counts")
count.matrix@x[count.matrix@x > 0] <- 1 #binarize counts
umap.coordinates <- Embeddings(object = KSAC, reduction = "umap.atac")

#update counts matrix
rownames <- rownames(count.matrix)
chrom <- sapply(rownames, function(x) strsplit(x, split = "-")[[1]][[1]])
loc_start <- sapply(rownames, function(x) strsplit(x, split = "-")[[1]][[2]])
loc_end <- sapply(rownames, function(x) strsplit(x, split = "-")[[1]][[3]])

unique(chrom)

pData <- data.frame(clusters)
fData <- data.frame(site_name = rownames(count.matrix), chromosome = chrom, bp1 = loc_start, bp2 = loc_end)
input_cds <- newCellDataSet(count.matrix, phenoData = new("AnnotatedDataFrame", data = pData),
                            featureData = new("AnnotatedDataFrame", data = fData),
                            expressionFamily = VGAM::binomialff())

cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = umap.coordinates, k = 50)
conns <- run_cicero(cicero_cds, hg38.chr.lengths) 

## Assign peaks to modules
ccan.assigns <- generate_ccans(conns)

links <- ConnectionsToLinks(conns = conns, ccans = ccan.assigns)
Links(KSAC) <- links




###





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

# snCv3 - Marker Genes and Comparison with Reference Data Sets (ED Fig 2)----------------------------
library(Seurat)
library(plyr)
library(dplyr)
library(Matrix)
library(ggplot2)
library("corrplot")

KBR <- LoadH5Seurat("Kidney_Healthy-Injury_Cell_Atlas_snCv3_Seurat_03282022.h5Seurat")


###Dotplot of marker gene expression (ED Fig 2a)
###Plot a subset of markers from the epi/int plots together
epi.order <- as.character(c(1:3,4:19,20,21,22,24,23,25,26,27,30,31,36,37,39,
                            28,29,32:35,38,40,41,44,45,
                            43,47,42,46,48,49,52,61,56,55,58,60,
                            50,51,53,54,57,59))
int.order <- as.character(c(62:65,67,70,
                            69,66,68,71:74,76:78,80:82,84:85, 
                            79,86,75,83,87,88:99))
KBR.sub <- subset(KBR, idents = c(epi.order,int.order))
Idents(object = KBR.sub) <- factor(Idents(object = KBR.sub), levels = c(epi.order,int.order))


epi.markers <- c("PTPRQ","WT1","NPHS2",                         #POD
                 "CDKN1C", "SPOCK2",                            #dPOD
                 "IGFBP2",                                      #dPOD
                 "CLDN1", "CFH","ALDH1A2",                      #PEC
                 "LRP2",                                        #PT
                 "SLC5A12","SLC22A6",                           #S1/S2
                 "PRODH2","SLC5A2","SLC22A8","SLC7A8",          #S1
                 "SLC34A1","SLC5A10",                           #S2                                  
                 "SLC5A11","SLC7A13",                           #S3
                 
                 "ITGB8","CDH6","HAVCR1","VCAM1",               #aPT
                 "CRYAB","TACSTD2","SLC44A5",                   #TL
                 "AQP1", "UNC5D",                               #DTL1
                 "ADGRL3","ID1",                                #DTL2
                 #DTL3
                 "AKR1B1","SH3GL3",                             #DTL3/ATL
                 "PROX1",                                       #ATL
                 
                 "CASR","SLC12A1","UMOD",                       #TAL
                 "PROM1",                                       #aTAL
                 "ITGB6","FGF13",                               #aTAL
                 "CLDN14", "KCTD16",                            #M-TAL
                 "ANK2",                                        #M-TAL
                 "ESRRB","EGF",                                 #TAL
                 "ENOX1","TMEM207",                             #C-TAL
                 "CLDN16",                                      #C-TAL
                 "NOS1","ROBO2",                                #MD
                 
                 "SLC12A3","TRPM6",                             #DCT
                 #DCT1
                 "SLC8A1","SCN2A","HSD11B2","CALB1",            #CNT
                 "TRPV5",                                       #DCT2
                 "SCNN1G","SCNN1B",                             #CNT                   
                 
                 "GATA3","AQP2","AQP3",                         #PC
                 "PDE10A",                                      #C-PC
                 "KCNK13",                                      #M-PC
                 "FXYD4",                                       #PC
                 "PHACTR1","SLC14A2",                           #IMCD
                 "TP63","GPX2",                                 #PapE
                 "ATP6V0D2",                                    #IC
                 "SLC4A1","SLC26A7",                            #IC-A
                 "KIT","AQP6","CALCA",                          #M-IC-A
                 "SLC4A9","INSRR")                              #IC-B


int.markers <- c("PECAM1","EMCN",                               #EC
                 "HECW2","PLAT","ITGA8",                        #EC-GC
                 "BTNL9","PALMD","AQP1","TM4SF1",               #EC-AEA-DVR
                 "SERPINE2",  #EC-AEA
                 "SLC14A1","ENPP2",                             #EC-DVR
                 "DNASE1L3","CEACAM1",                          #EC-AVR/PTC
                 "PITPNC1","SLCO2A1",                           #EC-PTC
                 "PLVAP","TLL1",                                #EC-AVR
                 "MMRN1","PROX1",                               #EC-LYM
                 
                 "PDGFRB",                                      #VSMC/P
                 
                 "ROBO1",                                       #MC/REN
                 "PIEZO2","POSTN",                              #MC
                 "REN","GRID2",                                 #REN
                 
                 "MYH11","RGS6","MCAM",                         #VSMC               
                 "RGS5","ADGRB3",                               #VSMC/P
                 
                 "COL1A1","COL1A2","C7","DCN",                  #FIB/MyoF
                 
                 "SYNPO2","PCDH7",                              #MyoF
                 "MEG3","LAMA2",                                #FIB/MyoF
                 "COL6A3","GLI2","COL5A1",                      #MyoF
                 "SLC24A3","CDH13","SMOC2","ANO3","RXFP1",      #MyoF
                 
                 "PDGFRA",                                      #FIB
                 
                 "SYT1",                                        #M-FIB
                 "PLCXD3","GREB1L",                             #M-FIB                          
                 "ADAMTSL1","FREM1",                            #M-FIB
                 "MGP","BGN",                                   #M-FIB
                 
                 
                 "ABI3BP","FLRT2", "FGF14","IGF1",              #aFIB
                 
                 "TAGLN","ACTA2",                               #dVSMC/P
                 "PTPRC",                                       #IMM
                 "BANK1","MS4A1",                               #B
                 "IGKC","MZB1",                                 #pB
                 
                 "THEMIS","IL7R","CD96","CD247",                #T
                 "GNLY","NKG7","GZMA",                          #NKT
                 
                 "KIT", "MS4A2",                                #Mast
                 
                 "MRC1","CD163","CD14",                         #MAC
                 "DIAPH3","CENPF","MKI67",                      #cycMAC
                 
                 "MSR1",
                 "ITGAX","HLA-DQA1","CSF2RA",                   #cDC
                 "FLT3","CLEC9A",                               #cDC
                 
                 "IL3RA","CLEC4C",                              #pDC
                 
                 "CTSS","FCN1","FCGR3A",                        #MON
                 "S100A9","S100A8","FCGR3B"                     #NC
                 
)

DotPlot(KBR.sub, features = unique(c(epi.markers,int.markers)),cols = c("lightgrey", "#225ea8"),dot.scale = 4) + RotatedAxis()






###Comparison with reference data sets
###Ransick et al. data (Ex Data Fig 2f) (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE129798)
mcm.meta <- readRDS("/media/Scratch_SSD_Voyager/Blue/KPMP_10X/KPMP_Ref_analysis_12-2020/McMahon_Mouse_2019/pmeta.rds")
mcm.data <- readRDS("/media/Scratch_SSD_Voyager/Blue/KPMP_10X/KPMP_Ref_analysis_12-2020/McMahon_Mouse_2019/expr_data.rds")
mcm.genes <- readRDS("/media/Scratch_SSD_Voyager/Blue/KPMP_10X/KPMP_Ref_analysis_12-2020/McMahon_Mouse_2019/gene_tbl.rds")
mcm.proj <- readRDS("/media/Scratch_SSD_Voyager/Blue/KPMP_10X/KPMP_Ref_analysis_12-2020/McMahon_Mouse_2019/kid_proj.rds")
mcm.otg <- readRDS("/media/Scratch_SSD_Voyager/Blue/KPMP_10X/KPMP_Ref_analysis_12-2020/McMahon_Mouse_2019/otg_tbl.rds")
mcm.otg <- data.frame(mcm.otg)

rownames(mcm.data) <- toupper(rownames(mcm.data))
mcm <- CreateSeuratObject(counts = mcm.data, project = "Mouse Kidney", min.cells = 3, min.features = 200)
mcm.meta <- na.omit(mcm.meta)
colnames(mcm)
mcm <- subset(mcm, cells = rownames(mcm.meta))
mcm <- NormalizeData(mcm, normalization.method = "LogNormalize", scale.factor = 10000)
mcm <- FindVariableFeatures(mcm, selection.method = "vst", nfeatures = 3000)
all.genes <- rownames(mcm)
mcm <- ScaleData(mcm, features = all.genes)

#Use only reference clusters
Idents(KBR) <- "state"
KBR.sub <- subset(KBR, idents = "reference")
Idents(KBR.sub) <- "clusters"
Idents(object = KBR.sub) <- factor(Idents(object = KBR.sub), levels = 1:119)
to.keep <- levels(Idents(KBR.sub))[levels(Idents(KBR.sub)) %in% 1:61]
KBR.sub <- subset(KBR.sub, idents = to.keep)
ave.KBR<-AverageExpression(KBR.sub, features = common.genes, slot = "scale.data")

#reorder Ransick et al. cell types
mcm@meta.data <- mcm.meta[rownames(mcm@meta.data),]
Idents(mcm) <- "Ontology_ID"
order <- c(1:8,"9A","9B",10:20,23,26,28:30,22,25,27,21,24,31,32)
Idents(object = mcm) <- factor(Idents(object = mcm), levels = order)

#Intersect all variable genes
common.genes <- intersect(VariableFeatures(KBR), VariableFeatures(mcm))

ave.KBR<-AverageExpression(KBR.sub, features = common.genes, slot = "scale.data")
ave.mcm<-AverageExpression(mcm, features = common.genes, slot = "scale.data")

ave.cor<-cor(cbind(ave.KBR$RNA,ave.mcm$RNA))
ave.cor<-ave.cor[1:26,27:59]

col2 <- colorRampPalette(c("#053061","#2166AC","#4393C3","#92C5DE","#D1E5F0","#FFFFFF","#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F"))
corrplot(ave.cor, col = col2(500),  method="color",order = "original", cl.lim=c(-1,1))



###Chen et al. scRNA-seq Data (Ex Data Fig 2c) (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE150338)
Idents(object = KBR) <- "clusters"
Idents(object = KBR) <- factor(Idents(object = KBR), levels = 1:119)
KBR.sub <- subset(KBR, idents = 1:100)

load("GSE150338_Seurat_IntegrateData.RData")
Knep.sc <- merge.data.integrated
DefaultAssay(Knep.sc) <- "RNA"
Idents(Knep.sc)
order <- c(10,4,5,13,12,8,0:3,6,9,7,11)
Idents(Knep.sc) <- factor(Idents(Knep.sc), levels = order)
Knep.sc <- ScaleData(Knep.sc)
DefaultAssay(Knep.sc) <- "integrated"

mouse.var <- VariableFeatures(Knep.sc)
DefaultAssay(Knep.sc) <- "RNA"

#Subset to nephron clusters
KBR.sub <- subset(KBR, idents = 26:59)

common.human.genes <- intersect(VariableFeatures(KBR.sub), toupper(mouse.var)) 
common.mouse.genes <- mouse.var[which(toupper(mouse.var) %in% common.human.genes)] 

#Use only reference clusters
ave.KBR <- AverageExpression(KBR.sub, features = common.human.genes, slot = "scale.data")
dim(ave.KBR$RNA)
ave.KBR$RNA <- ave.KBR$RNA[,c("40","41","44","45")]

ave.Kn <- AverageExpression(Knep.sc, features = common.mouse.genes, slot = "scale.data")
rownames(ave.Kn$integrated) <- toupper(rownames(ave.Kn$integrated))
ave.Kn$integrated <- ave.Kn$integrated[common.human.genes,c("0","1","2","3","6","9")]

ave.cor<-cor(cbind(ave.KBR$RNA,ave.Kn$integrated))
ave.cor<-ave.cor[1:4,5:10]
col2 <- colorRampPalette(c("#053061","#2166AC","#4393C3","#92C5DE","#D1E5F0","#FFFFFF","#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F"))
corrplot(ave.cor, col = col2(500),  method="color",order = "original", cl.lim=c(-1,1))




###Chen et al. Bulk Segment Data (Ex Data Fig 2d) (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE150338)
bulk <- read.delim("GSE150338_Mouse_TPM_Replicates.txt",sep="\t",header=TRUE, row.names = 1)
bulk <- bulk[!duplicated(bulk$Gene_symbol),]
rownames(bulk) <- bulk$Gene_symbol
bulk <- bulk[,3:110]

Kn.b <- CreateSeuratObject(counts = bulk, project = "Knepper", min.cells = 3, min.features = 100)
count.data <- GetAssayData(object = Kn.b[["RNA"]], slot = "counts")
Kn.b <- SetAssayData(object = Kn.b,  slot = "data", new.data = count.data, assay = "RNA")
all.genes <- rownames(Kn.b)
Kn.b <- ScaleData(Kn.b, features = all.genes)
Kn.b <- FindVariableFeatures(Kn.b, selection.method = "vst", nfeatures = 3000)
VariableFeatures(Kn.b)

common.human.genes <- intersect(VariableFeatures(KBR.sub), toupper(VariableFeatures(Kn.b))) 
common.mouse.genes <- VariableFeatures(Kn.b)[which(toupper(VariableFeatures(Kn.b)) %in% common.human.genes)] 

#Use only reference clusters
KBR.sub <- subset(KBR, idents = 1:61)

ave.KBR <- AverageExpression(KBR.sub, features = common.human.genes, slot = "scale.data")
order.1 <- c("4","5","6","7","20","21","22","24","31","36",
             "37","39","40","41","44","45","48","49","52","55","56","58","60","61")

ave.Kn <- AverageExpression(Kn.b, features = common.mouse.genes, slot = "scale.data")
rownames(ave.Kn$RNA) <- toupper(rownames(ave.Kn$RNA))
order.2 <- c("PTS1","PTS2","PTS3","DTL1","DTL2","DTL3","ATL","MTAL","CTAL","DCT","CNT","CCD",
             "OMCD","IMCD")

ave.cor<-cor(cbind(ave.KBR$RNA[rownames(ave.Kn$RNA),order.1],ave.Kn$RNA[,order.2]))
ave.cor<-ave.cor[1:24,25:38]

col2 <- colorRampPalette(c("#053061","#2166AC","#4393C3","#92C5DE","#D1E5F0","#FFFFFF","#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F"))
corrplot(ave.cor, col = col2(500),  method="color",order = "original", cl.lim=c(-1,1))



###Kirita et al. data (Ext Data Fig 2e) (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE139107)
meta <- read.table("GSE139107_MouseIRI.metadata.txt",sep="\t",header=TRUE,row.names=1)
H1 <- read.table("GSE139107_MouseIRI_control.dge.txt",sep="\t",header=TRUE,row.names=1)
H2 <- read.table("GSE139107_MouseIRI_4hours.dge.txt",sep="\t",header=TRUE,row.names=1)
H3 <- read.table("GSE139107_MouseIRI_12hours.dge.txt",sep="\t",header=TRUE,row.names=1)
H4 <- read.table("GSE139107_MouseIRI_2days.dge.txt",sep="\t",header=TRUE,row.names=1)
H5 <- read.table("GSE139107_MouseIRI_14days.dge.txt",sep="\t",header=TRUE,row.names=1)
H6 <- read.table("GSE139107_MouseIRI_6weeks.dge.txt",sep="\t",header=TRUE,row.names=1)
H1 <- as.sparse(H1)
H2 <- as.sparse(H2)
H3 <- as.sparse(H3)
H4 <- as.sparse(H4)
H5 <- as.sparse(H5)
H6 <- as.sparse(H6)

merge.sparse = function(listMatrixes) {
  # takes a list of sparse matrixes with different columns and adds them row wise
  
  allRownames <- sort(unique(unlist(lapply(listMatrixes,rownames))))
  allColnames <- as.character(unlist(lapply(listMatrixes,colnames)))
  for (currentMatrix in listMatrixes) {
    newRowLocations <- match(rownames(currentMatrix),allRownames)
    indexes <- which(currentMatrix>0, arr.ind = T)
    newRows <- newRowLocations[indexes[,1]]
    columns <- indexes[,2]
    newMatrix <- sparseMatrix(i=newRows,j=columns, x=currentMatrix@x,
                              dims=c(length(allRownames),max(columns)))
    if (!exists("matrixToReturn")) {
      matrixToReturn <- newMatrix
    }
    else {
      matrixToReturn <- cbind2(matrixToReturn,newMatrix)
    }
  }
  rownames(matrixToReturn) <- allRownames
  colnames(matrixToReturn) <- allColnames
  matrixToReturn  
}
H.list = mget(ls(pattern = "\\H"))
H.data <- merge.sparse(H.list)

rownames(H.data) <- toupper(rownames(H.data))
hki <- CreateSeuratObject(counts = H.data, project = "Mouse Kidney", min.cells = 3, min.features = 200,
                          meta.data = meta)
hki <- NormalizeData(hki, normalization.method = "LogNormalize", scale.factor = 10000)
hki <- FindVariableFeatures(hki, selection.method = "vst", nfeatures = 3000)
all.genes <- rownames(hki)
hki <- ScaleData(hki, features = all.genes)
Idents(hki) <- "celltype"
levels(Idents(hki))
Idents(hki) <- factor(Idents(hki), levels = c("Pod","PEC","PTS1","PTS2","PTS3","NewPT1","NewPT2",
                                              "DTL-ATL","MTAL","CTAL1","CTAL2","MD","DCT","DCT-CNT","CNT",
                                              "PC1","PC2","ICA","ICB","Uro","EC1","EC2","Per","Fib",         
                                              "Tcell","Mø"))
table(Idents(hki))
common.mouse <- intersect(rownames(KBR.sub), VariableFeatures(hki))  
common.human <- intersect(VariableFeatures(KBR.sub)[1:3000], rownames(hki)) 
common.genes <- unique(c(common.mouse, common.human))

#Correlate at the cluster level
Idents(object = KBR.sub) <- "clusters"
Idents(object = KBR.sub) <- factor(Idents(object = KBR.sub), levels = 1:100)

common.mouse <- intersect(rownames(KBR.sub), VariableFeatures(hki))  
common.human <- intersect(VariableFeatures(KBR.sub)[1:3000], rownames(hki)) 
common.genes <- unique(c(common.mouse, common.human))

ave.KBR<-AverageExpression(KBR.sub, features = common.genes, slot = "scale.data")
ave.hki<-AverageExpression(hki, features = common.genes, slot = "scale.data")

ave.cor<-cor(cbind(ave.KBR$RNA,ave.hki$RNA))
ave.cor<-ave.cor[1:12,101:107]

col2 <- colorRampPalette(c("#053061","#2166AC","#4393C3","#92C5DE","#D1E5F0","#FFFFFF","#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F"))
corrplot(ave.cor, col = col2(500),  method="color",order = "original",cl.lim=c(-1,1))




###Kuppe et al. (Ext Data Fig 2g) (https://zenodo.org/record/4059315#.YApbyehKiHs)
Idents(object = KBR) <- "clusters"
Idents(object = KBR) <- factor(Idents(object = KBR), levels = 1:119)
KBR.sub <- subset(KBR, idents = 1:100)

##CD10+ data 
colData <- read.table("Human_CD10plus/PTmap_UMI_counts_colData.txt",sep=",",header=TRUE)
rowData <- read.table("Human_CD10plus/PTmap_UMI_counts_rowData.txt",sep=",",header=TRUE)
Data <- readMM("Human_CD10plus/PTmap_UMI_counts.mtx")
colnames(Data) <- rownames(colData)
rownames(Data) <- rowData$Gene.Symbol

cd10p <- CreateSeuratObject(counts = Data, project = "PT", min.cells = 3, min.features = 200,
                            meta.data = colData)

##CD10- data 
colData <- read.table("CD10negative/kidneyMap_UMI_counts_colData.txt",sep=",",header=TRUE)
rowData <- read.table("CD10negative/kidneyMap_UMI_counts_rowData.txt",sep=",",header=TRUE)
Data <- readMM("CD10negative/kidneyMap_UMI_counts.mtx")
colnames(Data) <- rownames(colData)
rownames(Data) <- rowData$Gene.Symbol

cd10n <- CreateSeuratObject(counts = Data, min.cells = 3, min.features = 200,
                            meta.data = colData)
cd10n$Annotation.Level.3 <- as.character(cd10n$Annotation.Level.3)
cd10n$Annotation.Level.3[cd10n$Annotation.Level.3 == "Fibroblast 6"] <- "Fibroblast 3"
cd10n$Annotation.Level.3[cd10n$Annotation.Level.3 == "Fibroblast 2"] <- "Fibroblast 1"
cd10n$Annotation.Level.3[cd10n$Annotation.Level.3 == "Fibroblast 4"] <- "Fibroblast 2"
table(cd10n$Annotation.Level.3)

#PDGFRb data 
colData <- read.table("Human_PDGFRb/pdgfrbMap_UMI_counts_colData.txt",sep=",",header=TRUE)
rowData <- read.table("Human_PDGFRb/pdgfrbMap_UMI_counts_rowData.txt",sep=",",header=TRUE)
Data <- readMM("Human_PDGFRb/pdgfrbMap_UMI_counts.mtx")
colnames(Data) <- rownames(colData)
rownames(Data) <- rowData$Gene.Symbol

PDGFRb <- CreateSeuratObject(counts = Data, min.cells = 3, min.features = 200,
                             meta.data = colData)
PDGFRb$Annotation.Level.3 <- as.character(PDGFRb$Annotation.Level.3)
PDGFRb$Annotation.Level.3[PDGFRb$Annotation.Level.3 == "Fibroblast 6"] <- "Fibroblast 3"
PDGFRb$Annotation.Level.3[PDGFRb$Annotation.Level.3 == "Fibroblasts 4a"] <- "Fibroblast 2a"
PDGFRb$Annotation.Level.3[PDGFRb$Annotation.Level.3 == "Fibroblast 2"] <- "Fibroblast 1"
PDGFRb$Annotation.Level.3[PDGFRb$Annotation.Level.3 == "Fibroblasts 4b"] <- "Fibroblast 2b"
table(PDGFRb$Annotation.Level.3)


Kr.Seurat <- merge(x = cd10p, y = c(cd10n, PDGFRb))
Idents(Kr.Seurat) <- "Annotation.Level.3"
Kr.Seurat <- NormalizeData(Kr.Seurat, normalization.method = "LogNormalize", scale.factor = 10000)
Kr.Seurat <- FindVariableFeatures(Kr.Seurat, selection.method = "vst", nfeatures = 3000)
all.genes <- rownames(Kr.Seurat)
Kr.Seurat <- ScaleData(Kr.Seurat, features = all.genes)
Idents(Kr.Seurat) <- "Annotation.Level.2"
order <- c("Podocytes","PT","Proximal Tubule","Injured Proximal Tubule",
           "Descending Thin Limb","Thick Ascending Limb","Injured Tubule",
           "Distal Convoluted Tubule","Connecting Tubule","Collecting Duct Principal Cells",
           "Prinicipal Cells","Intercalated Cells","Uroethlial Cells",
           "Glomerular Capillaries","Arterioles","Arterioral Endothelium","Vasa Recta",                   
           "Venular Endothelium","Lymph Endothelium","Injured Endothelium",            
           "Mesangial Cells","Pericytes","Vascular Smooth Muscle Cells",
           "Smooth Muscle Cells","Schwann Cells","Myofibroblast","Myofibroblasts",                   
           "Fibroblast","Fibroblasts",
           "B Cells","Plasma Cells","T Cells","Natural Killer Cells",
           "Macrophages","Dendritic Cells",                                
           "Monocytes","Basophils"                      
)
Idents(Kr.Seurat) <- factor(Idents(Kr.Seurat), levels = order)

#Subset to only the stromal cell types (cortical only)
KBR.st <- subset(KBR.sub, idents = c(71:81,86:87))
Idents(KBR.st) <- factor(Idents(KBR.st), levels = c(71:81,86:87))
Idents(Kr.Seurat) <- "Annotation.Level.3"
levels(Idents(Kr.Seurat))
order <- c("Mesangial Cells 1","Vascular Smooth Muscle Cells","Vascular Smooth Muscle Cells 2",
           "Pericytes","Pericytes 1","Pericytes 2",
           "Myofibroblasts 1","Myofibroblast 1a","Myofibroblast 1b","Myofibroblasts 2a",
           "Myofibroblasts 2b","Myofibroblasts 3a","Myofibroblasts 3b",                                        
           "Fibroblast 1","Fibroblast 2","Fibroblast 2a","Fibroblast 2b","Fibroblast 3")                                                  

Kr.Seurat.st <- subset(Kr.Seurat, idents = order)
Kr.Seurat.st <- FindVariableFeatures(Kr.Seurat.st, selection.method = "vst", nfeatures = 3000)
Idents(Kr.Seurat.st) <- factor(Idents(Kr.Seurat.st), levels = order)

common.genes <- intersect(VariableFeatures(KBR.st), VariableFeatures(Kr.Seurat.st)) 

ave.KBR <- AverageExpression(KBR.st, features = common.genes, slot = "scale.data")
ave.Kr <- AverageExpression(Kr.Seurat.st, features = common.genes, slot = "scale.data")

ave.cor<-cor(cbind(ave.KBR$RNA,ave.Kr$RNA))
ave.cor<-ave.cor[1:13,14:31]

col2 <- colorRampPalette(c("#053061","#2166AC","#4393C3","#92C5DE","#D1E5F0","#FFFFFF","#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F"))
corrplot(ave.cor, col = col2(500),  method="color",order = "original",cl.lim=c(-1,1))




###Immune cell types (Ext Data Fig 2h-i)

#Prepare immune reference data
library(celldex)
library(SingleR)

ref.se <- MonacoImmuneData()

library(scater)
library(Seurat)
library(SeuratDisk)
library(dplyr)
library(Matrix)
library(patchwork)
library(ggplot2)
library("corrplot")

setwd("/media/Scratch_SSD_Voyager/Blue/KPMP_10X/KPMP_Ref_analysis_12-2020")
load("Kidney_KPMP-Biopsy_Ref_10X-R_092020_Seurat_D.rda")
Idents(object = KBR) <- "clusters"
Idents(object = KBR) <- factor(Idents(object = KBR), levels = 1:119)
KBR.sub <- subset(KBR, idents = 1:100)
KBR.imm <- subset(KBR.sub, idents = 88:99)
Idents(KBR.imm) <- factor(Idents(KBR.imm), levels = c(88:99))



##Compare with human sorted immune populations (PBMC) (Ext Data Fig 2i)
logcounts <- assays(ref.se)$logcounts
MID <- CreateSeuratObject(counts = logcounts, min.cells = 0, min.features = 0)
MID@meta.data$label.main <- colData(ref.se)$label.main
MID@meta.data$label.fine <- colData(ref.se)$label.fine
MID@meta.data$label.ont <- colData(ref.se)$label.ont

Idents(MID) <- "label.main"
order <- c("B cells","CD4+ T cells","CD8+ T cells","T cells","NK cells",             
           "Monocytes","Dendritic cells","Neutrophils","Progenitors","Basophils")
Idents(MID) <- factor(Idents(MID), levels = order)

MID <- SetAssayData(
  object = MID,
  slot = "data",
  new.data = logcounts,
  assay = "RNA"
)
all.genes <- rownames(MID)
MID <- ScaleData(MID, features = all.genes)

#Identify and remove poorly aligned subtypes using fine labels
Idents(MID) <- "label.fine"
all.genes <- rownames(MID)

#Use Immune Genes (https://www.immport.org/shared/genelists)
immune.genes <- read.delim("GeneList.txt")
common.genes <- intersect(intersect(rownames(KBR.imm), immune.genes$Symbol),rownames(MID))

#Use only reference clusters
ave.KBR<-AverageExpression(KBR.imm, features = common.genes, slot = "scale.data")
ave.MID<-AverageExpression(MID, features = common.genes, slot = "scale.data")

ave.cor<-cor(cbind(ave.KBR$RNA,ave.MID$RNA))
ave.cor<-ave.cor[1:12,13:41]
ave.cor<-ave.cor[,colSums(ave.cor) > 0.2]

#Remove cell types with poor correlation
MID <- subset(MID, idents = colnames(ave.cor))

Idents(MID) <- "label.main"
order <- c("B cells","CD4+ T cells","CD8+ T cells","T cells","NK cells",             
           "Monocytes","Dendritic cells","Neutrophils","Progenitors","Basophils")
Idents(MID) <- factor(Idents(MID), levels = order)

#Use Immune Genes
common.genes <- intersect(intersect(rownames(KBR.imm), immune.genes$Symbol),rownames(MID))

ave.KBR<-AverageExpression(KBR.imm, features = common.genes, slot = "scale.data")
ave.MID<-AverageExpression(MID, features = common.genes, slot = "scale.data")
ave.cor<-cor(cbind(ave.KBR$RNA,ave.MID$RNA))
ave.cor<-ave.cor[1:12,13:19]
col2 <- colorRampPalette(c("#053061","#2166AC","#4393C3","#92C5DE","#D1E5F0","#FFFFFF","#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F"))
corrplot(ave.cor, col = col2(500),  method="color",order = "original", cl.lim=c(-1,1))



##Compare with mouse immune populations (Ext Data Fig 2h)
library(celldex)
library(SingleR)

ref <- ImmGenData()

logcounts <- assays(ref)$logcounts
rownames(logcounts) <- toupper(rownames(logcounts))

IGD <- CreateSeuratObject(counts = logcounts, min.cells = 0, min.features = 0)
IGD@meta.data$label.main <- colData(ref)$label.main
IGD@meta.data$label.fine <- colData(ref)$label.fine
IGD@meta.data$label.ont <- colData(ref)$label.ont

Idents(IGD) <- "label.main"
order <- c("Epithelial cells","Endothelial cells","Stromal cells","Fibroblasts","Microglia",
           "B cells","B cells, pro","T cells","NKT","NK cells","Mast cells",
           "Macrophages","Monocytes","DC","Neutrophils","Eosinophils","Basophils",             
           "ILC","Tgd","Stem cells")
Idents(IGD) <- factor(Idents(IGD), levels = order)

IGD <- SetAssayData(
  object = IGD,
  slot = "data",
  new.data = logcounts,
  assay = "RNA"
)
all.genes <- rownames(IGD)
IGD <- ScaleData(IGD, features = all.genes)

IGD <- subset(IGD, idents = c("B cells","B cells, pro","T cells","NKT","NK cells","Mast cells",
                              "Macrophages","Monocytes","DC","Neutrophils"))

#Identify and remove poorly aligned subtypes using fine labels
Idents(IGD) <- "label.fine"
all.genes <- rownames(IGD)

#Use Immune Genes (https://www.immport.org/shared/genelists)
immune.genes <- read.delim("GeneList.txt")
common.genes <- intersect(intersect(rownames(KBR.imm), immune.genes$Symbol),rownames(IGD))

ave.KBR<-AverageExpression(KBR.imm, features = common.genes, slot = "scale.data")
ave.IGD<-AverageExpression(IGD, features = common.genes, slot = "scale.data")

ave.cor<-cor(cbind(ave.KBR$RNA,ave.IGD$RNA))
ave.cor<-ave.cor[1:12,13:191]
ave.cor<-ave.cor[,colSums(ave.cor) > 0.2]

#Remove cell types with poor correlation
IGD <- subset(IGD, idents = colnames(ave.cor))

Idents(IGD) <- "label.main"
order <- c("B cells","B cells, pro","T cells","NKT","NK cells","Mast cells",
           "Macrophages","Monocytes","DC","Neutrophils")
Idents(IGD) <- factor(Idents(IGD), levels = order)

common.genes <- intersect(intersect(rownames(KBR.imm), immune.genes$Symbol),rownames(IGD))
ave.KBR<-AverageExpression(KBR.imm, features = common.genes, slot = "scale.data")
ave.IGD<-AverageExpression(IGD, features = common.genes, slot = "scale.data")

ave.cor<-cor(cbind(ave.KBR$RNA,ave.IGD$RNA))
ave.cor<-ave.cor[1:12,13:20]

col2 <- colorRampPalette(c("#053061","#2166AC","#4393C3","#92C5DE","#D1E5F0","#FFFFFF","#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F"))
corrplot(ave.cor, col = col2(500),  method="color",order = "original",
         hclust.method = "ward.D", cl.lim=c(-1,1))




###Predicted Lake et al. 2019 (Nat. Comm.) snDrop-seq Annotations (ED Fig 2b)
##Load snDrop data
#https://drive.google.com/drive/folders/1p-0PNmLDnaTTvz39fr3hA1WKo0rvCWzc?usp=sharing
load("Adult_Normal_Kidney_snDrop_Lake2019_NCOMM_Seuratv3.Robj") 
#meta <- read.table("Kidney_snDrop_Meta.txt",sep="\t",header=TRUE,row.names=1)
#ank.3@meta.data <- meta[rownames(ank.3@meta.data),]
#ank.3$assay <- "snDrop"

##Use variable features as integration anchors
ref.p2 #pagoda2 object generated when clustering snDrop data (not provided)
ref.var.info <- ref.p2$misc$varinfo; ref.var.info <- ref.var.info[order(ref.var.info$lp, decreasing = F),];
ref.var.genes <- rownames(ref.var.info)[1:2000]
query.p2 #pagoda2 object generated when clustering snCv3 data (not provided)
q.var.info <- query.p2$misc$varinfo; q.var.info <- q.var.info[order(q.var.info$lp, decreasing = F),];
query.var.genes <- rownames(q.var.info)[1:2000]

genes.use <- unique(c(ref.var.genes, query.var.genes))
genes.use <- intersect(genes.use, intersect(rownames(ank.3),rownames(KBR))) #2603

NK.anchors <- FindTransferAnchors(reference = ank.3, query = KBR, features = genes.use, 
                                  reduction = "cca")
celltype.predictions <- TransferData(anchorset = NK.anchors, refdata = ank.3$Subclass, 
                                     weight.reduction = "cca")

###

# snCv3 - Degenerative State Gene Set and Score (Supplementary Table 12) -----------------------------------------
library(Seurat)
library(SeuratDisk)
library(dplyr)
library(Matrix)
library(ggplot2)
library(forcats)
library(cowplot)
library(patchwork)

KBR <- LoadH5Seurat("Kidney_Healthy-Injury_Cell_Atlas_snCv3_Seurat_03282022.h5Seurat")
Idents(object = KBR) <- "clusters"
Idents(object = KBR) <- factor(Idents(object = KBR), levels = 1:119)
KBR.sub <- subset(KBR, idents = 1:100)

load("color_factors.robj")


#Load Function
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


#Prepare Seurat object
Idents(KBR.sub) <- "subclass.l3"
Idents(KBR.sub) <- factor(Idents(KBR.sub), levels = c(
  "POD","dPOD","PEC","PT-S1","PT-S2","PT-S3","aPT","cycPT","dPT","dPT/DTL",
  "DTL2","DTL1","DTL3","dDTL3","ATL","dATL","aTAL1","aTAL2","M-TAL","dM-TAL",
  "C-TAL","dC-TAL","MD","DCT1","DCT2","dDCT","cycDCT","CNT","CNT-PC","dCNT",
  "cycCNT","CCD-PC","OMCD-PC","dOMCD-PC","IMCD","dIMCD","CCD-IC-A","CNT-IC-A",
  "dC-IC-A","OMCD-IC-A","tPC-IC","IC-B","PapE","EC-GC","EC-AEA","EC-DVR","EC-PTC",
  "dEC-PTC","EC-AVR","dEC","cycEC","EC-LYM","MC","REN","VSMC","VSMC/P","dVSMC",
  "MYOF","cycMYOF","FIB","M-FIB","dM-FIB","aFIB","dFIB","B","PL","T","NKT","MAST",
  "MAC-M2","cycMNP","MDC","cDC","pDC","ncMON","N","SC/NEU"
))
table(Idents(KBR.sub))





###Conserved Degen State Markers - POD
KBR.sub.pod <- subset(KBR.sub, idents = c("POD","dPOD"))

celltype <- levels(Idents(KBR.sub.pod))

cons.POD.DEGs <- getconservedDEGs(KBR.sub.pod, celltype)



###Conserved Degen State Markers - PT S1-S3
KBR.sub.pts13 <- subset(KBR.sub, idents = c("PT-S1","PT-S2","PT-S3",
                                            "dPT","dPT/DTL"))
KBR.sub.pts13 <- RenameIdents(object = KBR.sub.pts13, 
                              "PT-S1" = "PT", 
                              "PT-S2" = "PT",
                              "PT-S3" = "PT",
                              "dPT/DTL" = "dPT")
KBR.sub.pts13$subclass.l3 <- Idents(KBR.sub.pts13)

celltype <- levels(Idents(KBR.sub.pts13))
cons.pts13.DEGs <- getconservedDEGs(KBR.sub.pts13, celltype)



###Conserved Degen State Markers - DTL3 - ATL
KBR.sub.dtl3 <- subset(KBR.sub, idents = c("DTL3", "dDTL3", "ATL", "dATL"))
KBR.sub.dtl3 <- RenameIdents(object = KBR.sub.dtl3, 
                             "DTL3" = "DTL-ATL",
                             "dDTL3" = "dDTL-ATL",
                             "ATL" = "DTL-ATL",
                             "dATL" = "dDTL-ATL")

celltype <- levels(Idents(KBR.sub.dtl3))
cons.dtl3.DEGs <- getconservedDEGs(KBR.sub.dtl3, celltype)



###Conserved Degen State Markers - TAL
KBR.sub.tal <- subset(KBR.sub, idents = c("M-TAL","dM-TAL","C-TAL","dC-TAL","MD"))
KBR.sub.tal <- RenameIdents(object = KBR.sub.tal, 
                            "M-TAL" = "TAL",
                            "dM-TAL" = "dTAL",
                            "C-TAL" = "TAL",
                            "dC-TAL" = "dTAL",
                            "MD" = "TAL")
celltype <- levels(Idents(KBR.sub.tal))
cons.tal.DEGs <- getconservedDEGs(KBR.sub.tal, celltype)



###Conserved Degen State Markers - DCT
KBR.sub.dct <- subset(KBR.sub, idents = c("DCT1","DCT2","dDCT"))
KBR.sub.dct <- RenameIdents(object = KBR.sub.dct, 
                            "DCT1" = "DCT",
                            "DCT2" = "DCT")
KBR.sub.dct$subclass.l3 <- Idents(KBR.sub.dct)

celltype <- levels(Idents(KBR.sub.dct))
cons.dct.DEGs <- getconservedDEGs(KBR.sub.dct, celltype)




###Conserved Degen State Markers - CNT
KBR.sub.cnt <- subset(KBR.sub, idents = c("CNT","CNT-PC","dCNT"))
KBR.sub.cnt <- RenameIdents(object = KBR.sub.cnt, 
                            "CNT-PC" = "CNT")
celltype <- levels(Idents(KBR.sub.cnt))
cons.cnt.DEGs <- getconservedDEGs(KBR.sub.cnt, celltype)




###Conserved Degen State Markers - PC  (not including IMCD)
KBR.sub.pc <- subset(KBR.sub, idents = c("CCD-PC","OMCD-PC","dOMCD-PC"))
KBR.sub.pc <- RenameIdents(object = KBR.sub.pc, 
                           "CCD-PC" = "PC",
                           "OMCD-PC" = "PC",
                           "dOMCD-PC" = "dPC")
celltype <- levels(Idents(KBR.sub.pc))
cons.pc.DEGs <- getconservedDEGs(KBR.sub.pc, celltype)




###Conserved Degen State Markers - IC-A
KBR.sub.ic <- subset(KBR.sub, idents = c("CCD-IC-A","CNT-IC-A","dC-IC-A","OMCD-IC-A"))
KBR.sub.ic <- RenameIdents(object = KBR.sub.ic, 
                           "CCD-IC-A" = "IC-A",
                           "CNT-IC-A" = "IC-A",
                           "dC-IC-A" = "dIC-A",
                           "OMCD-IC-A" = "IC-A")
celltype <- levels(Idents(KBR.sub.ic))
cons.ic.DEGs <- getconservedDEGs(KBR.sub.ic, celltype)



###Conserved Degen State Markers - EC
KBR.sub.ec <- subset(KBR.sub, idents = c("EC-GC","EC-AEA","EC-DVR","EC-PTC","dEC-PTC",
                                         "EC-AVR","dEC"))
KBR.sub.ec <- RenameIdents(object = KBR.sub.ec, 
                           "EC-GC" = "EC",
                           "EC-AEA" = "EC",
                           "EC-DVR" = "EC",
                           "EC-PTC" = "EC",
                           "dEC-PTC" = "dEC",
                           "EC-AVR" = "EC",
                           "dEC" = "dEC")
celltype <- levels(Idents(KBR.sub.ec))
cons.ec.DEGs <- getconservedDEGs(KBR.sub.ec, celltype)



###Conserved Degen State Markers - VSMC/P
KBR.sub.smc <- subset(KBR.sub, idents = c("MC","REN","VSMC","VSMC/P","dVSMC"))
KBR.sub.smc <- RenameIdents(object = KBR.sub.smc, 
                            "MC" = "VSMC",
                            "REN" = "VSMC",
                            "VSMC" = "VSMC",
                            "VSMC/P" = "VSMC",
                            "dVSMC" = "dVSMC")
celltype <- levels(Idents(KBR.sub.smc))
cons.smc.DEGs <- getconservedDEGs(KBR.sub.smc, celltype)



###Conserved Degen State Markers - FIB - All
KBR.sub.fib <- subset(KBR.sub, idents = c("FIB","M-FIB","dM-FIB","dFIB"))
KBR.sub.fib <- RenameIdents(object = KBR.sub.fib, 
                            "M-FIB" = "FIB",
                            "dM-FIB" = "dFIB")
celltype <- levels(Idents(KBR.sub.fib))
cons.fib.DEGs <- getconservedDEGs(KBR.sub.fib, celltype)




###Combine markers for conserved gene set
DEG.list = list(cons.POD.DEGs,
                cons.pts13.DEGs,
                cons.dtl3.DEGs,
                cons.tal.DEGs,
                cons.dct.DEGs,
                cons.cnt.DEGs,
                cons.pc.DEGs,
                cons.ic.DEGs,
                cons.ec.DEGs,
                cons.smc.DEGs,
                cons.fib.DEGs)
celltype = 1:11

cons.DEGs <- do.call(rbind, lapply(celltype, function(ct) {
  print(paste("Running for cell type:", ct))
  
  cl.mark <- DEG.list[[ct]]
  cl.mark <- cl.mark[cl.mark$AKI_p_val < 0.05 & 
                       cl.mark$CKD_p_val < 0.05 & 
                       cl.mark$Ref_p_val < 0.05,]
  return(cl.mark)
  
}))




###Add in Degen score
degen.states <- c("dPOD", "dPT", "dDTL-ATL", "dTAL", "dDCT",  "dCNT", "dPC", "dIC-A",
                  "dEC", "dVSMC", "dFIB")
degen.gene.counts <- data.frame(cons.DEGs[cons.DEGs$cluster %in% degen.states,] %>%
                                  group_by(gene) %>%
                                  tally())
degen.gene.counts <- degen.gene.counts[order(-degen.gene.counts$n),]
degen.genes <- degen.gene.counts[degen.gene.counts$n >= 4, ]$gene #154
count.data = GetAssayData(KBR, slot = "data")
count.data = as.matrix(count.data[rownames(count.data) %in% toupper(degen.genes) ,])
combined_matrix = rbind(count.data,average=apply(count.data,2,mean))
cor_matrix = cor(t(combined_matrix))
cor_vector = cor_matrix[,dim(cor_matrix)[1]]
count.data_restricted = count.data[rownames(count.data) %in% names(cor_vector[cor_vector >= 0.1]),]
aaa = apply(count.data_restricted,2,mean)

KBR$degen.score <- aaa[rownames(KBR@meta.data)]



###Combine markers for expanded gene set
DEG.list = list(cons.POD.DEGs,
                cons.pts13.DEGs,
                cons.dtl3.DEGs,
                cons.tal.DEGs,
                cons.dct.DEGs,
                cons.cnt.DEGs,
                cons.pc.DEGs,
                cons.ic.DEGs,
                cons.ec.DEGs,
                cons.smc.DEGs,
                cons.fib.DEGs)
celltype = 1:11

cons.DEGs <- do.call(rbind, lapply(celltype, function(ct) {
  print(paste("Running for cell type:", ct))
  
  cl.mark <- DEG.list[[ct]]
  cl.mark <- subset(cl.mark, AKI_avg_log2FC >= 0.6 | CKD_avg_log2FC >= 0.6 | Ref_avg_log2FC >= 0.6)
  cl.mark <- subset(cl.mark, AKI_p_val_adj < 0.05 | CKD_p_val_adj < 0.05 | Ref_p_val_adj < 0.05)
  cl.mark <- cl.mark[cl.mark$max_pval < 0.05,]
  return(cl.mark)
  
}))




# snCv3 - Adaptive Epithelial State Gene Set and Score (Supplementary Table 13) -----------------------------------------
library(Seurat)
library(SeuratDisk)
library(dplyr)
library(Matrix)
library(ggplot2)
library(viridis)
library(swne)

KBR <- LoadH5Seurat("Kidney_Healthy-Injury_Cell_Atlas_snCv3_Seurat_03282022.h5Seurat")
Idents(object = KBR) <- "clusters"
Idents(object = KBR) <- factor(Idents(object = KBR), levels = 1:119)
KBR.sub <- subset(KBR, idents = 1:100)


###Conserved aEpi Markers

#Load Function
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

#Prepare Seurat object
Idents(KBR.sub) <- "subclass.l3"
Idents(KBR.sub) <- factor(Idents(KBR.sub), levels = c(
  "POD","dPOD","PEC","PT-S1","PT-S2","PT-S3","aPT","cycPT","dPT","dPT/DTL",
  "DTL2","DTL1","DTL3","dDTL3","ATL","dATL","aTAL1","aTAL2","M-TAL","dM-TAL",
  "C-TAL","dC-TAL","MD","DCT1","DCT2","dDCT","cycDCT","CNT","CNT-PC","dCNT",
  "cycCNT","CCD-PC","OMCD-PC","dOMCD-PC","IMCD","dIMCD","CCD-IC-A","CNT-IC-A",
  "dC-IC-A","OMCD-IC-A","tPC-IC","IC-B","PapE","EC-GC","EC-AEA","EC-DVR","EC-PTC",
  "dEC-PTC","EC-AVR","dEC","cycEC","EC-LYM","MC","REN","VSMC","VSMC/P","dVSMC",
  "MYOF","cycMYOF","FIB","M-FIB","dM-FIB","aFIB","dFIB","B","PL","T","NKT","MAST",
  "MAC-M2","cycMNP","MDC","cDC","pDC","ncMON","N","SC/NEU"))
table(Idents(KBR.sub))



###Conserved Markers of aPT
KBR.sub.pt <- subset(KBR.sub, idents = c("PT-S1","PT-S2","PT-S3","aPT"))
KBR.sub.pt <- RenameIdents(object = KBR.sub.pt, 
                           "PT-S1" = "PT", 
                           "PT-S2" = "PT",
                           "PT-S3" = "PT")
celltype <- levels(Idents(KBR.sub.pt))
cons.pt.as.DEGs <- getconservedDEGs(KBR.sub.pt, celltype)



###Conserved Markers of aTAL
KBR.sub.tal <- subset(KBR.sub, idents = c("aTAL1","M-TAL","C-TAL"))
KBR.sub.tal <- RenameIdents(object = KBR.sub.tal, 
                            "aTAL1" = "aTAL", 
                            "M-TAL" = "TAL",
                            "C-TAL" = "TAL")
table(Idents(KBR.sub.tal))
Idents(KBR.sub.tal) <- factor(Idents(KBR.sub.tal), levels = c("TAL","aTAL"))
celltype <- levels(Idents(KBR.sub.tal))
cons.tal.as.DEGs <- getconservedDEGs(KBR.sub.tal, celltype)





###Conserved markers across aPT and aTAL
###Combine markers
DEG.list = list(cons.pt.as.DEGs,
                cons.tal.as.DEGs)
celltype = 1:2

cons.DEGs <- do.call(rbind, lapply(celltype, function(ct) {
  print(paste("Running for cell type:", ct))
  
  cl.mark <- DEG.list[[ct]]
  cl.mark <- cl.mark[cl.mark$AKI_p_val < 0.05 & 
                       cl.mark$CKD_p_val < 0.05 & 
                       cl.mark$Ref_p_val < 0.05,]
  top10 <- cl.mark %>% 
    mutate(avg_log2FC = (AKI_avg_log2FC + CKD_avg_log2FC + Ref_avg_log2FC) /3) %>% 
    group_by(cluster) %>% 
    top_n(n = 10, 
          wt = avg_log2FC)
  return(cl.mark)
  
}))




###Select conserved markers that are more specific to the aEpi state 
KBR.sub <- RenameIdents(KBR.sub, "aPT" = "aTub")
KBR.sub <- RenameIdents(KBR.sub, "aTAL1" = "aTub")
KBR.sub <- RenameIdents(KBR.sub, "aTAL2" = "aTub")

aTub.markers <- FindMarkers(KBR.sub, ident.1 = "aTub", features = common.adap.genes,
                            only.pos = TRUE)



###Altered state score
adap.state.genes <- rownames(aTub.markers) #82 genes
count.data = GetAssayData(KBR, slot = "data")
count.data = as.matrix(count.data[rownames(count.data) %in% toupper(adap.state.genes) ,])
combined_matrix = rbind(count.data,average=apply(count.data,2,mean))
cor_matrix = cor(t(combined_matrix))
cor_vector = cor_matrix[,dim(cor_matrix)[1]]
count.data_restricted = count.data[rownames(count.data) %in% names(cor_vector[cor_vector >= 0.1]),]
aaa = apply(count.data_restricted,2,mean)
KBR$aEpi.score <- aaa[rownames(KBR@meta.data)]




###Expanded adaptive -epi state gene set
DEG.list = list(cons.pt.as.DEGs,
                cons.tal.as.DEGs)
celltype = 1:2

cons.DEGs <- do.call(rbind, lapply(celltype, function(ct) {
  print(paste("Running for cell type:", ct))
  
  cl.mark <- DEG.list[[ct]]
  cl.mark <- subset(cl.mark, AKI_avg_log2FC >= 0.6 | CKD_avg_log2FC >= 0.6 | Ref_avg_log2FC >= 0.6)
  cl.mark <- subset(cl.mark, AKI_p_val_adj < 0.05 | CKD_p_val_adj < 0.05 | Ref_p_val_adj < 0.05)
  cl.mark <- cl.mark[cl.mark$max_pval < 0.05,]
  return(cl.mark)
  
}))





###Gene Set Enrichment Analysis (GSEA) (Ext Data Fig. 8c-d)
GSEA = function(gene_list, GO_file, pval) {
  set.seed(54321)
  library(dplyr)
  library(gage)
  library(fgsea)
  library(data.table)
  if ( any( duplicated(names(gene_list)) )  ) {
    warning("Duplicates in gene names")
    gene_list = gene_list[!duplicated(names(gene_list))]
  }
  if  ( !all( order(gene_list, decreasing = TRUE) == 1:length(gene_list)) ){
    warning("Gene list not sorted")
    gene_list = sort(gene_list, decreasing = TRUE)
  }
  myGO = fgsea::gmtPathways(GO_file)
  
  fgRes <- fgsea::fgsea(pathways = myGO, 
                        stats = gene_list,
                        minSize=10,
                        maxSize=500,
                        eps=0) %>% 
    as.data.frame() #%>% 
  #dplyr::filter(padj < !!pval)
  print(dim(fgRes))
  
  ## Filter FGSEA by using gage results. Must be significant and in same direction to keep 
  gaRes = gage::gage(gene_list, gsets=myGO, same.dir=TRUE, set.size =c(10,500))
  
  ups = as.data.frame(gaRes$greater) %>% 
    tibble::rownames_to_column("Pathway") %>% 
    dplyr::filter(!is.na(p.geomean) ) %>%
    dplyr::select("Pathway")
  
  downs = as.data.frame(gaRes$less) %>% 
    tibble::rownames_to_column("Pathway") %>% 
    dplyr::filter(!is.na(p.geomean) ) %>%
    dplyr::select("Pathway")
  
  print(dim(rbind(ups,downs)))
  keepups = as.data.table(fgRes[fgRes$NES > 0 & !is.na(match(fgRes$pathway, ups$Pathway)), ])
  keepdowns = as.data.table(fgRes[fgRes$NES < 0 & !is.na(match(fgRes$pathway, downs$Pathway)), ])
  
  ### Collapse redundant pathways
  Up = fgsea::collapsePathways(keepups[order(pval)], pathways = myGO, stats = gene_list,  nperm = 500, pval.threshold = 0.05)
  Down = fgsea::collapsePathways(keepdowns[order(pval)], myGO, gene_list,  nperm = 500, pval.threshold = 0.05) 
  
  fgRes = fgRes[ !is.na(match(fgRes$pathway, 
                              c( Up$mainPathways, Down$mainPathways))), ] %>% 
    arrange(desc(NES))
  fgRes$pathway = stringr::str_replace(fgRes$pathway, "GO_" , "")
  
  fgRes$Enrichment = ifelse(fgRes$NES > 0, "Up-regulated", "Down-regulated")
  filtRes = rbind(head(fgRes, n = 10),
                  tail(fgRes, n = 10 ))
  g = ggplot(filtRes, aes(reorder(pathway, NES), NES)) +
    geom_segment( aes(reorder(pathway, NES), xend=pathway, y=0, yend=NES)) +
    geom_point( size=5, aes( fill = Enrichment),
                shape=21, stroke=2) +
    scale_fill_manual(values = c("Down-regulated" = "dodgerblue",
                                 "Up-regulated" = "firebrick") ) +
    coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score",
         title="GSEA - Biological Processes") + 
    theme_minimal()
  
  output = list("Results" = fgRes, "Plot" = g)
  return(output)
} #(from https://bioinformaticsbreakdown.com/how-to-gsea/)
GO_file <- "c5.all.v7.2.symbols.gmt"  #GO set (http://www.gsea-msigdb.org/gsea/downloads.jsp#msigdb)


#aPT
S4table <- cons.pt.as.DEGs
S4table$avg_log2FC <- (as.numeric(S4table$AKI_avg_log2FC) + as.numeric(S4table$CKD_avg_log2FC) + as.numeric(S4table$Ref_avg_log2FC)) /3
S4table[S4table$cluster == "aPT",]$avg_log2FC <- as.numeric(S4table[S4table$cluster == "aPT",]$avg_log2FC)
S4table[S4table$cluster == "PT",]$avg_log2FC <- (-1 * as.numeric(S4table[S4table$cluster == "PT",]$avg_log2FC))
gene_list = as.numeric(S4table$avg_log2FC)
names(gene_list) = as.character(S4table$gene)
gene_list = sort(gene_list, decreasing = TRUE)
gene_list = gene_list[!duplicated(names(gene_list))]
res = GSEA(gene_list, GO_file, pval = 0.05)
res$Plot


#aTAL
S4table <- cons.tal.as.DEGs
S4table$avg_log2FC <- (as.numeric(S4table$AKI_avg_log2FC) + as.numeric(S4table$CKD_avg_log2FC) + as.numeric(S4table$Ref_avg_log2FC)) /3
S4table[S4table$cluster == "aTAL",]$avg_log2FC <- as.numeric(S4table[S4table$cluster == "aTAL",]$avg_log2FC)
S4table[S4table$cluster == "TAL",]$avg_log2FC <- (-1 * as.numeric(S4table[S4table$cluster == "TAL",]$avg_log2FC))
gene_list = as.numeric(S4table$avg_log2FC)
names(gene_list) = as.character(S4table$gene)
gene_list = sort(gene_list, decreasing = TRUE)
gene_list = gene_list[!duplicated(names(gene_list))]
res = GSEA(gene_list, GO_file, pval = 0.05)
res$Plot






# snCv3 - Adaptive Stromal State Gene Set and Score (Supplementary Table 14) -------------------------------------
library(Seurat)
library(SeuratDisk)
library(dplyr)
library(Matrix)
library(ggplot2)
library(viridis)
library(swne)

KBR <- LoadH5Seurat("Kidney_Healthy-Injury_Cell_Atlas_snCv3_Seurat_03282022.h5Seurat")
Idents(object = KBR) <- "clusters"
Idents(object = KBR) <- factor(Idents(object = KBR), levels = 1:119)
KBR.sub <- subset(KBR, idents = 1:100)



###Conserved Adaptive Stromal Markers
#Load Function
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


#Prepare Seurat object
Idents(KBR.sub) <- "subclass.l3"
Idents(KBR.sub) <- factor(Idents(KBR.sub), levels = c(
  "POD","dPOD","PEC","PT-S1","PT-S2","PT-S3","aPT","cycPT","dPT","dPT/DTL",
  "DTL2","DTL1","DTL3","dDTL3","ATL","dATL","aTAL1","aTAL2","M-TAL","dM-TAL",
  "C-TAL","dC-TAL","MD","DCT1","DCT2","dDCT","cycDCT","CNT","CNT-PC","dCNT",
  "cycCNT","CCD-PC","OMCD-PC","dOMCD-PC","IMCD","dIMCD","CCD-IC-A","CNT-IC-A",
  "dC-IC-A","OMCD-IC-A","tPC-IC","IC-B","PapE","EC-GC","EC-AEA","EC-DVR","EC-PTC",
  "dEC-PTC","EC-AVR","dEC","cycEC","EC-LYM","MC","REN","VSMC","VSMC/P","dVSMC",
  "MYOF","cycMYOF","FIB","M-FIB","dM-FIB","aFIB","dFIB","B","PL","T","NKT","MAST",
  "MAC-M2","cycMNP","MDC","cDC","pDC","ncMON","N","SC/NEU"))




###Conserved Markers of aFIB
KBR.sub.fib <- subset(KBR.sub, idents = c("FIB","aFIB"))
KBR.sub.fib$subclass.l3 <- Idents(KBR.sub.fib)
celltype <- levels(Idents(KBR.sub.fib))
cons.fib.as.DEGs <- getconservedDEGs(KBR.sub.fib, celltype)




###Conserved Markers of MYOF
KBR.sub.fib <- subset(KBR.sub, idents = c("MYOF","FIB"))
Idents(KBR.sub.fib) <- factor(Idents(KBR.sub.fib), levels = c("FIB", "MYOF"))
KBR.sub.fib$subclass.l3 <- Idents(KBR.sub.fib)
celltype <- levels(Idents(KBR.sub.fib))
cons.myof.as.DEGs <- getconservedDEGs(KBR.sub.fib, celltype)





###Conserved Markers for aFIB and MYOF combined
KBR.sub.fib <- subset(KBR.sub, idents = c("MYOF","FIB","aFIB"))
KBR.sub.fib <- RenameIdents(KBR.sub.fib, 'MYOF' = "aFIB")
Idents(KBR.sub.fib) <- factor(Idents(KBR.sub.fib), levels = c("FIB", "aFIB"))
KBR.sub.fib$subclass.l3 <- Idents(KBR.sub.fib)
celltype <- levels(Idents(KBR.sub.fib))
cons.aStr.DEGs <- getconservedDEGs(KBR.sub.fib, celltype)


###Conserved markers across aStr
###Combine markers
DEG.list = list(cons.aStr.DEGs)
celltype = 1

cons.DEGs <- do.call(rbind, lapply(celltype, function(ct) {
  print(paste("Running for cell type:", ct))
  
  cl.mark <- DEG.list[[ct]]
  cl.mark <- cl.mark[cl.mark$AKI_p_val < 0.05 & 
                       cl.mark$CKD_p_val < 0.05 & 
                       cl.mark$Ref_p_val < 0.05,]
  
  return(cl.mark)
  
}))



###Select conserved markers that are more specific to the altered state 
KBR.sub <- RenameIdents(KBR.sub, "MYOF" = "aFIB")
aFIB.markers <- FindMarkers(KBR.sub, ident.1 = "aFIB", features = aStr.genes,
                            only.pos = TRUE)




###Altered state score
alt.state.genes <- rownames(aFIB.markers[aFIB.markers$avg_log2FC > 0.6,]) #42 genes
count.data = GetAssayData(KBR, slot = "data")
count.data = as.matrix(count.data[rownames(count.data) %in% toupper(alt.state.genes) ,])
combined_matrix = rbind(count.data,average=apply(count.data,2,mean))
cor_matrix = cor(t(combined_matrix))
cor_vector = cor_matrix[,dim(cor_matrix)[1]]
count.data_restricted = count.data[rownames(count.data) %in% names(cor_vector[cor_vector >= 0.1]),]
aaa = apply(count.data_restricted,2,mean)
KBR$aStr.state.score <- aaa[rownames(KBR@meta.data)]




###Expanded adaptive -str state gene set
DEG.list = list(cons.fib.as.DEGs, cons.myof.as.DEGs)
celltype = 1:2

cons.DEGs <- do.call(rbind, lapply(celltype, function(ct) {
  print(paste("Running for cell type:", ct))
  
  cl.mark <- DEG.list[[ct]]
  cl.mark <- subset(cl.mark, AKI_avg_log2FC >= 0.6 | CKD_avg_log2FC >= 0.6 | Ref_avg_log2FC >= 0.6)
  cl.mark <- subset(cl.mark, AKI_p_val_adj < 0.05 | CKD_p_val_adj < 0.05 | Ref_p_val_adj < 0.05)
  cl.mark <- cl.mark[cl.mark$max_pval < 0.05,]
  return(cl.mark)
  
}))








# snCv3 - Cycling State Gene Set and Score (Supplementary Table 15) ------------------------------------------
library(Seurat)
library(SeuratDisk)
library(dplyr)
library(Matrix)
library(ggplot2)
library(forcats)
library(cowplot)
library(patchwork)

KBR <- LoadH5Seurat("Kidney_Healthy-Injury_Cell_Atlas_snCv3_Seurat_03282022.h5Seurat")
Idents(object = KBR) <- "clusters"
Idents(object = KBR) <- factor(Idents(object = KBR), levels = 1:119)
KBR.sub <- subset(KBR, idents = 1:100)


#Load Function
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


#Prepare Seurat object
Idents(KBR.sub) <- "subclass.l3"
Idents(KBR.sub) <- factor(Idents(KBR.sub), levels = c(
  "POD","dPOD","PEC","PT-S1","PT-S2","PT-S3","aPT","cycPT","dPT","dPT/DTL",
  "DTL2","DTL1","DTL3","dDTL3","ATL","dATL","aTAL1","aTAL2","M-TAL","dM-TAL",
  "C-TAL","dC-TAL","MD","DCT1","DCT2","dDCT","cycDCT","CNT","CNT-PC","dCNT",
  "cycCNT","CCD-PC","OMCD-PC","dOMCD-PC","IMCD","dIMCD","CCD-IC-A","CNT-IC-A",
  "dC-IC-A","OMCD-IC-A","tPC-IC","IC-B","PapE","EC-GC","EC-AEA","EC-DVR","EC-PTC",
  "dEC-PTC","EC-AVR","dEC","cycEC","EC-LYM","MC","REN","VSMC","VSMC/P","dVSMC",
  "MYOF","cycMYOF","FIB","M-FIB","dM-FIB","aFIB","dFIB","B","PL","T","NKT","MAST",
  "MAC-M2","cycMNP","MDC","cDC","pDC","ncMON","N","SC/NEU")
)



###Cycling Markers - PT
KBR.sub.pts13 <- subset(KBR.sub, idents = c("PT-S1","PT-S2","PT-S3","cycPT"))
KBR.sub.pts13 <- RenameIdents(object = KBR.sub.pts13, 
                              "PT-S1" = "PT", 
                              "PT-S2" = "PT",
                              "PT-S3" = "PT")
pts13.DEGs <- FindAllMarkers(KBR.sub.pts13, min.pct = 0.25, only.pos = TRUE, 
                             max.cells.per.ident = 300)



###Cycling Markers - DCT
KBR.sub.dct <- subset(KBR.sub, idents = c("DCT1","DCT2","cycDCT"))
KBR.sub.dct <- RenameIdents(object = KBR.sub.dct, 
                            "DCT1" = "DCT",
                            "DCT2" = "DCT")
dct.DEGs <- FindAllMarkers(KBR.sub.dct, min.pct = 0.25, only.pos = TRUE, 
                           max.cells.per.ident = 300)




###Cycling Markers - CNT
KBR.sub.cnt <- subset(KBR.sub, idents = c("CNT","CNT-PC","cycCNT"))
KBR.sub.cnt <- RenameIdents(object = KBR.sub.cnt, 
                            "CNT-PC" = "CNT")
cnt.DEGs <- FindAllMarkers(KBR.sub.cnt, min.pct = 0.25, only.pos = TRUE, 
                           max.cells.per.ident = 300)



###Cycling Markers - EC
KBR.sub.ec <- subset(KBR.sub, idents = c("EC-GC","EC-AEA","EC-DVR","EC-PTC",
                                         "EC-AVR","cycEC"))
KBR.sub.ec <- RenameIdents(object = KBR.sub.ec, 
                           "EC-GC" = "EC",
                           "EC-AEA" = "EC",
                           "EC-DVR" = "EC",
                           "EC-PTC" = "EC",
                           "EC-AVR" = "EC")
ec.DEGs <- FindAllMarkers(KBR.sub.ec, min.pct = 0.25, only.pos = TRUE, 
                          max.cells.per.ident = 300)




###Cycling markers - MYOF
KBR.sub.fib <- subset(KBR.sub, idents = c("MYOF","cycMYOF"))

fib.DEGs <- FindAllMarkers(KBR.sub.fib, min.pct = 0.25, only.pos = TRUE, 
                           max.cells.per.ident = 300)




###Combine markers for conserved gene set
DEG.list = list(pts13.DEGs,
                dct.DEGs,
                cnt.DEGs,
                ec.DEGs,
                fib.DEGs)
celltype = 1:5

cons.DEGs <- do.call(rbind, lapply(celltype, function(ct) {
  print(paste("Running for cell type:", ct))
  
  cl.mark <- DEG.list[[ct]]
  cl.mark <- cl.mark[cl.mark$p_val_adj < 0.05 & cl.mark$avg_log2FC > 0.6,]
  return(cl.mark)
  
}))




###Cycling state score
cyc.state.genes <- cyc.gene.counts[cyc.gene.counts$n == 5, ]$gene #69
count.data = GetAssayData(KBR, slot = "data")
count.data = as.matrix(count.data[rownames(count.data) %in% toupper(cyc.state.genes) ,])
combined_matrix = rbind(count.data,average=apply(count.data,2,mean))
cor_matrix = cor(t(combined_matrix))
cor_vector = cor_matrix[,dim(cor_matrix)[1]]
count.data_restricted = count.data[rownames(count.data) %in% names(cor_vector[cor_vector >= 0.1]),]
aaa = apply(count.data_restricted,2,mean)
KBR$cyc.score <- aaa[rownames(KBR@meta.data)]




###Add in cell cycle score using Seurat
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
KBR <- CellCycleScoring(KBR, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)





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

# scCv3 - Degenerative State Gene Set and Score (Supplementary Table 12) --------------------------------------------------
library(Seurat)
library(SeuratDisk)
library(dplyr)
library(Matrix)
library(ggplot2)
library(forcats)
library(cowplot)
library(patchwork)

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


#Load Function
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


###Conserved Degenerative Markers - POD
sc.KBR.pod <- subset(sc.KBR, idents = c("POD","dPOD"))
celltype <- levels(Idents(sc.KBR.pod))
cons.POD.DEGs <- getconservedDEGs(sc.KBR.pod, celltype)



###Conserved Degen Markers - PT S1-S3
sc.KBR.pts13 <- subset(sc.KBR, idents = c("PT-S1/2","PT-S3",
                                          "dPT"))
sc.KBR.pts13 <- RenameIdents(object = sc.KBR.pts13, 
                             "PT-S1/2" = "PT", 
                             "PT-S3" = "PT"
)
celltype <- levels(Idents(sc.KBR.pts13))
cons.pts13.DEGs <- getconservedDEGs(sc.KBR.pts13, celltype)



###Conserved Degen Markers - TAL
sc.KBR.tal <- subset(sc.KBR, idents = c("M-TAL","dM-TAL","C-TAL","dC-TAL","MD"))
sc.KBR.tal <- RenameIdents(object = sc.KBR.tal, 
                           "M-TAL" = "TAL",
                           "dM-TAL" = "dTAL",
                           "C-TAL" = "TAL",
                           "dC-TAL" = "dTAL",
                           "MD" = "TAL")
celltype <- levels(Idents(sc.KBR.tal))
cons.tal.DEGs <- getconservedDEGs(sc.KBR.tal, celltype)



###Conserved Degen Markers - DCT
sc.KBR.dct <- subset(sc.KBR, idents = c("DCT1","DCT2","dDCT"))
sc.KBR.dct <- RenameIdents(object = sc.KBR.dct, 
                           "DCT1" = "DCT",
                           "DCT2" = "DCT")
celltype <- levels(Idents(sc.KBR.dct))
cons.dct.DEGs <- getconservedDEGs(sc.KBR.dct, celltype)



###Conserved Degen Markers - CNT
sc.KBR.cnt <- subset(sc.KBR, idents = c("CNT","CNT-PC","dCNT"))
sc.KBR.cnt <- RenameIdents(object = sc.KBR.cnt, 
                           "CNT-PC" = "CNT")
celltype <- levels(Idents(sc.KBR.cnt))
cons.cnt.DEGs <- getconservedDEGs(sc.KBR.cnt, celltype)



###Conserved Degen Markers - PC  (not including IMCD)
sc.KBR.pc <- subset(sc.KBR, idents = c("CCD-PC","OMCD-PC","dOMCD-PC"))
sc.KBR.pc <- RenameIdents(object = sc.KBR.pc, 
                          "CCD-PC" = "PC",
                          "OMCD-PC" = "PC",
                          "dOMCD-PC" = "dPC")
celltype <- levels(Idents(sc.KBR.pc))
cons.pc.DEGs <- getconservedDEGs(sc.KBR.pc, celltype)



###Conserved Degen Markers - IC-A
sc.KBR.ic <- subset(sc.KBR, idents = c("CCD-IC-A","CNT-IC-A","dC-IC-A","OMCD-IC-A"))
sc.KBR.ic <- RenameIdents(object = sc.KBR.ic, 
                          "CCD-IC-A" = "IC-A",
                          "CNT-IC-A" = "IC-A",
                          "dC-IC-A" = "dIC-A",
                          "OMCD-IC-A" = "IC-A")
celltype <- levels(Idents(sc.KBR.ic))
cons.ic.DEGs <- getconservedDEGs(sc.KBR.ic, celltype)



###Conserved Degen Markers - VSMC/P
sc.KBR.smc <- subset(sc.KBR, idents = c("MC","VSMC","VSMC/P","dVSMC"))
sc.KBR.smc <- RenameIdents(object = sc.KBR.smc, 
                           "MC" = "VSMC",
                           "VSMC" = "VSMC",
                           "VSMC/P" = "VSMC",
                           "dVSMC" = "dVSMC")
celltype <- levels(Idents(sc.KBR.smc))
cons.smc.DEGs <- getconservedDEGs(sc.KBR.smc, celltype)



###Conserved Degen Markers - FIB 
sc.KBR.fib <- subset(sc.KBR, idents = c("FIB","M-FIB","dFIB"))
sc.KBR.fib <- RenameIdents(object = sc.KBR.fib, 
                           "M-FIB" = "FIB")
celltype <- levels(Idents(sc.KBR.fib))
cons.fib.DEGs <- getconservedDEGs(sc.KBR.fib, celltype)



###Combine markers
DEG.list = list(cons.POD.DEGs,
                cons.pts13.DEGs,
                cons.tal.DEGs,
                cons.dct.DEGs,
                cons.cnt.DEGs,
                cons.pc.DEGs,
                cons.ic.DEGs,
                cons.smc.DEGs,
                cons.fib.DEGs)
celltype = 1:9

cons.DEGs <- do.call(rbind, lapply(celltype, function(ct) {
  print(paste("Running for cell type:", ct))
  
  cl.mark <- DEG.list[[ct]]
  cl.mark <- cl.mark[cl.mark$AKI_p_val < 0.05 & 
                       cl.mark$CKD_p_val < 0.05 & 
                       cl.mark$Ref_p_val < 0.05,]
  top10 <- cl.mark %>% 
    mutate(avg_log2FC = (AKI_avg_log2FC + CKD_avg_log2FC + Ref_avg_log2FC) /3) %>% 
    group_by(cluster) %>% 
    top_n(n = 10, 
          wt = avg_log2FC)
  return(cl.mark)
  
}))




###Add in degen score
degen.states <- c("dPOD", "dPT", "dTAL", "dDCT",  "dCNT", "dPC", "dIC-A",
                  "dVSMC", "dFIB")
degen.gene.counts <- data.frame(cons.DEGs[cons.DEGs$cluster %in% degen.states,] %>%
                                  group_by(gene) %>%
                                  tally())
degen.gene.counts <- degen.gene.counts[order(-degen.gene.counts$n),]
degen.mark <- degen.gene.counts[degen.gene.counts$n >= 4, ]$gene
count.data = GetAssayData(sc.KBR, slot = "data")
count.data = as.matrix(count.data[rownames(count.data) %in% toupper(degen.mark) ,])
combined_matrix = rbind(count.data,average=apply(count.data,2,mean))
cor_matrix = cor(t(combined_matrix))
cor_vector = cor_matrix[,dim(cor_matrix)[1]]
count.data_restricted = count.data[rownames(count.data) %in% names(cor_vector[cor_vector >= 0.1]),]
aaa = apply(count.data_restricted,2,mean)
sc.KBR$degen.score <- aaa[rownames(sc.KBR@meta.data)]




###Combine markers for expanded gene set
DEG.list = list(cons.POD.DEGs,
                cons.pts13.DEGs,
                cons.tal.DEGs,
                cons.dct.DEGs,
                cons.cnt.DEGs,
                cons.pc.DEGs,
                cons.ic.DEGs,
                cons.smc.DEGs,
                cons.fib.DEGs)
celltype = 1:9

cons.DEGs <- do.call(rbind, lapply(celltype, function(ct) {
  print(paste("Running for cell type:", ct))
  
  cl.mark <- DEG.list[[ct]]
  cl.mark <- subset(cl.mark, AKI_avg_log2FC >= 0.6 | CKD_avg_log2FC >= 0.6 | Ref_avg_log2FC >= 0.6)
  cl.mark <- subset(cl.mark, AKI_p_val_adj < 0.05 | CKD_p_val_adj < 0.05 | Ref_p_val_adj < 0.05)
  cl.mark <- cl.mark[cl.mark$max_pval < 0.05,]
  return(cl.mark)
  
}))



# scCv3 - Adaptive Epithelial State Gene Set and Score (Supplementary Table 13) ---------------------------------------------------------------
library(Seurat)
library(SeuratDisk)
library(dplyr)
library(Matrix)
library(ggplot2)
library(forcats)
library(cowplot)
library(patchwork)

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


#Load Function
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




###Conserved Markers of aPT
KBR.sub.pt <- subset(sc.KBR, idents = c("PT-S1/2","PT-S3","aPT"))
KBR.sub.pt <- RenameIdents(object = KBR.sub.pt, 
                           "PT-S1/2" = "PT",
                           "PT-S3" = "PT")
table(Idents(KBR.sub.pt))
celltype <- levels(Idents(KBR.sub.pt))
cons.pt.as.DEGs <- getconservedDEGs(KBR.sub.pt, celltype)



###Conserved Markers of aTAL
KBR.sub.tal <- subset(sc.KBR, idents = c("aTAL1","M-TAL","C-TAL"))
KBR.sub.tal <- RenameIdents(object = KBR.sub.tal, 
                            "aTAL1" = "aTAL", 
                            "M-TAL" = "TAL",
                            "C-TAL" = "TAL")
Idents(KBR.sub.tal) <- factor(Idents(KBR.sub.tal), levels = c("TAL","aTAL"))
celltype <- levels(Idents(KBR.sub.tal))
cons.tal.as.DEGs <- getconservedDEGs(KBR.sub.tal, celltype)




###Conserved aEpi markers
###Combine markers
DEG.list = list(cons.pt.as.DEGs,
                cons.tal.as.DEGs)
celltype = 1:2

cons.DEGs <- do.call(rbind, lapply(celltype, function(ct) {
  print(paste("Running for cell type:", ct))
  
  cl.mark <- DEG.list[[ct]]
  cl.mark <- cl.mark[cl.mark$AKI_p_val < 0.05 & 
                       cl.mark$CKD_p_val < 0.05 & 
                       cl.mark$Ref_p_val < 0.05,]
  top10 <- cl.mark %>% 
    mutate(avg_log2FC = (AKI_avg_log2FC + CKD_avg_log2FC + Ref_avg_log2FC) /3) %>% 
    group_by(cluster) %>% 
    top_n(n = 10, 
          wt = avg_log2FC)
  return(cl.mark)
  
}))

alt.states <- c("aPT","aTAL")
alt.gene.counts <- data.frame(cons.DEGs[cons.DEGs$cluster %in% alt.states,] %>%
                                group_by(gene) %>%
                                tally())
alt.gene.counts <- alt.gene.counts[order(-alt.gene.counts$n),]
common.alt.genes <- alt.gene.counts[alt.gene.counts$n == 2, ]$gene



###Select conserved markers that are more specific to the aEpi state 
sc.KBR <- RenameIdents(sc.KBR, "aPT" = "aTub")
sc.KBR <- RenameIdents(sc.KBR, "aTAL1" = "aTub")
sc.KBR <- RenameIdents(sc.KBR, "aTAL2" = "aTub")
aTub.markers <- FindMarkers(sc.KBR, ident.1 = "aTub", features = common.alt.genes,
                            only.pos = TRUE)




###Adaptive state score
adap.state.genes <- rownames(aTub.markers) #63 genes
adap.state.genes <- adap.state.genes[!adap.state.genes %in% "MALAT1"]
count.data = GetAssayData(sc.KBR, slot = "data")
count.data = as.matrix(count.data[rownames(count.data) %in% toupper(adap.state.genes) ,])
combined_matrix = rbind(count.data,average=apply(count.data,2,mean))
cor_matrix = cor(t(combined_matrix))
cor_vector = cor_matrix[,dim(cor_matrix)[1]]
count.data_restricted = count.data[rownames(count.data) %in% names(cor_vector[cor_vector >= 0.1]),]
aaa = apply(count.data_restricted,2,mean)
sc.KBR$aEpi.score <- aaa[rownames(sc.KBR@meta.data)]




###Expanded adaptive -epi state gene set
###Conserved markers across aPT and aTAL
###Combine markers

DEG.list = list(cons.pt.as.DEGs,
                cons.tal.as.DEGs)
celltype = 1:2

cons.DEGs <- do.call(rbind, lapply(celltype, function(ct) {
  print(paste("Running for cell type:", ct))
  
  cl.mark <- DEG.list[[ct]]
  cl.mark <- subset(cl.mark, AKI_avg_log2FC >= 0.6 | CKD_avg_log2FC >= 0.6 | Ref_avg_log2FC >= 0.6)
  cl.mark <- subset(cl.mark, AKI_p_val_adj < 0.05 | CKD_p_val_adj < 0.05 | Ref_p_val_adj < 0.05)
  cl.mark <- cl.mark[cl.mark$max_pval < 0.05,]
  return(cl.mark)
  
}))




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




# scCv3 - Cycling State Gene Set and Score (Supplementary Table 15) --------------------------------------------------
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


#Load Function
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


###Cycling Markers - PT
KBR.sub.pts13 <- subset(sc.KBR, idents = c("PT-S1/2","PT-S3","cycPT"))
KBR.sub.pts13 <- RenameIdents(object = KBR.sub.pts13, 
                              "PT-S1/2" = "PT",
                              "PT-S3" = "PT")

pts13.DEGs <- FindAllMarkers(KBR.sub.pts13, min.pct = 0.25, only.pos = TRUE, 
                             max.cells.per.ident = 300)




###Cycling Markers - EC
KBR.sub.ec <- subset(sc.KBR, idents = c("EC-GC","EC-AEA/DVR","EC-PTC",
                                        "EC-AVR","cycEC"))
KBR.sub.ec <- RenameIdents(object = KBR.sub.ec, 
                           "EC-GC" = "EC",
                           "EC-AEA/DVR" = "EC",
                           "EC-PTC" = "EC",
                           "EC-AVR" = "EC")
ec.DEGs <- FindAllMarkers(KBR.sub.ec, min.pct = 0.25, only.pos = TRUE, 
                          max.cells.per.ident = 300)




###Combine markers for conserved gene set
DEG.list = list(pts13.DEGs,
                ec.DEGs)
celltype = 1:2

cons.DEGs <- do.call(rbind, lapply(celltype, function(ct) {
  print(paste("Running for cell type:", ct))
  
  cl.mark <- DEG.list[[ct]]
  cl.mark <- cl.mark[cl.mark$p_val_adj < 0.05 & cl.mark$avg_log2FC > 0.6,]
  return(cl.mark)
  
}))

cyc.states <- c("cycPT","cycEC")
cyc.gene.counts <- data.frame(cons.DEGs[cons.DEGs$cluster %in% cyc.states,] %>%
                                group_by(gene) %>%
                                tally())
cyc.gene.counts <- cyc.gene.counts[order(-cyc.gene.counts$n),]





###Cycling state score
cyc.state.genes <- cyc.gene.counts[cyc.gene.counts$n == 2, ]$gene #101
count.data = GetAssayData(sc.KBR, slot = "data")
count.data = as.matrix(count.data[rownames(count.data) %in% toupper(cyc.state.genes) ,])
combined_matrix = rbind(count.data,average=apply(count.data,2,mean))
cor_matrix = cor(t(combined_matrix))
cor_vector = cor_matrix[,dim(cor_matrix)[1]]
count.data_restricted = count.data[rownames(count.data) %in% names(cor_vector[cor_vector >= 0.1]),]
aaa = apply(count.data_restricted,2,mean)
sc.KBR$cycling.score <- aaa[rownames(sc.KBR@meta.data)]







# snCv3/scCv3/SNARE - Altered State Signatures (ED Fig 8) -------------------------------------------------------------
library(Seurat)
library(SeuratDisk)
library(dplyr)
library(Matrix)
library(ggplot2)

load("sc-sn_int_color_factors.robj")
load("color_factors.robj")

###Prepare seurat objects
KBR <- LoadH5Seurat("Kidney_Healthy-Injury_Cell_Atlas_snCv3_Seurat_03282022.h5Seurat")
Idents(object = KBR) <- "clusters"
Idents(object = KBR) <- factor(Idents(object = KBR), levels = 1:119)
KBR <- subset(KBR, idents = 1:100)
KBR$id <- "ref"
sc.KBR <- LoadH5Seurat("Kidney_Healthy-Injury_Cell_Atlas_scCv3_Seurat_03282022.h5Seurat")
colnames(sc.KBR@meta.data)[colnames(sc.KBR@meta.data) == "state"] <- "state.l1"

#Combined snCv3/scCv3 Object
refquery <- merge(KBR, y = sc.KBR)
umap.coordinates <- rbind(Embeddings(object = KBR, reduction = "umap"),
                          Embeddings(object = sc.KBR, reduction = "ref.umap"))
refquery[["umap"]] <- CreateDimReducObject(embeddings = umap.coordinates, key = "umap_", assay = DefaultAssay(refquery))
refquery <- NormalizeData(refquery, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(refquery)
refquery <- ScaleData(refquery, features = all.genes, split.by = "id", assay = "RNA")
select.markers <- intersect(VariableFeatures(KBR), rownames(refquery))
VariableFeatures(refquery) <- select.markers
Idents(refquery) <- "subclass.l3"


###Plot Figures
###snCv3 - Degenerative State Features (Ext Data Fig 8a)
KBR.sub <- subset(KBR, idents = 1:100)
Idents(object = KBR.sub) <- "state.l2"
KBR.sub.inj <- subset(KBR.sub, idents = c("reference","degenerative"))
Idents(object = KBR.sub.inj) <- "state.l1"
Idents(object = KBR.sub.inj) <- factor(Idents(object = KBR.sub.inj), levels = c("reference","altered"))
KBR.sub.inj$state.l1 <- Idents(object = KBR.sub.inj)

Idents(KBR.sub.inj) <- "subclass.l1"
Idents(object = KBR.sub.inj) <- factor(Idents(object = KBR.sub.inj), 
                                       levels = c("POD","PEC","PT","DTL","ATL","TAL",
                                                  "DCT","CNT","PC","IC","PapE","EC","VSM/P",
                                                  "FIB","IMM","NEU"))
KBR.sub.inj <- subset(KBR.sub.inj, idents = c("POD","PT","DTL","ATL","TAL",
                                              "DCT","CNT","PC","IC","EC","VSM/P",
                                              "FIB"))

VlnPlot(object = KBR.sub.inj, features = c("IGFBP7","CLU","CST3","degen.score","nFeature_RNA","percent.er","percent.mt"), 
        ncol = 1, pt.size = -1, split.by = "state.l1", cols = as.character(state.l1.cols))



###Matrisome Score Plots (Ext Data Fig 8g)
#core matrisome score
matrisome_set = read.table("https://raw.githubusercontent.com/mahmoudibrahim/KidneyMap/master/assets/public/ecm_genes_human.txt", sep = "\t", header = TRUE)
matrisome = matrisome_set$Gene.Symbol[which(matrisome_set$Division %in% c("Core matrisome"))]
count.data = GetAssayData(KBR, slot = "data")
count.data = as.matrix(count.data[rownames(count.data) %in% toupper(matrisome) ,])
combined_matrix = rbind(count.data,average=apply(count.data,2,mean))
cor_matrix = cor(t(combined_matrix))
cor_vector = cor_matrix[,dim(cor_matrix)[1]]
count.data_restricted = count.data[rownames(count.data) %in% names(cor_vector[cor_vector >= 0.1]),]
aaa = apply(count.data_restricted,2,mean)

KBR$matrisome.score <- aaa[rownames(KBR@meta.data)]

#collagen score
matrisome_set = read.table("https://raw.githubusercontent.com/mahmoudibrahim/KidneyMap/master/assets/public/ecm_genes_human.txt", sep = "\t", header = TRUE)
matrisome = matrisome_set$Gene.Symbol[which(matrisome_set$Category %in% c("Collagens"))]
count.data = GetAssayData(KBR, slot = "data")
count.data = as.matrix(count.data[rownames(count.data) %in% toupper(matrisome) ,])
combined_matrix = rbind(count.data,average=apply(count.data,2,mean))
cor_matrix = cor(t(combined_matrix))
cor_vector = cor_matrix[,dim(cor_matrix)[1]]
count.data_restricted = count.data[rownames(count.data) %in% names(cor_vector[cor_vector >= 0.1]),]
aaa = apply(count.data_restricted,2,mean)

KBR$collagen.score <- aaa[rownames(KBR@meta.data)]

#glycoprotein score
matrisome_set = read.table("https://raw.githubusercontent.com/mahmoudibrahim/KidneyMap/master/assets/public/ecm_genes_human.txt", sep = "\t", header = TRUE)
matrisome = matrisome_set$Gene.Symbol[which(matrisome_set$Category %in% c("ECM Glycoproteins"))]
count.data = GetAssayData(KBR, slot = "data")
count.data = as.matrix(count.data[rownames(count.data) %in% toupper(matrisome) ,])
combined_matrix = rbind(count.data,average=apply(count.data,2,mean))
cor_matrix = cor(t(combined_matrix))
cor_vector = cor_matrix[,dim(cor_matrix)[1]]
count.data_restricted = count.data[rownames(count.data) %in% names(cor_vector[cor_vector >= 0.1]),]
aaa = apply(count.data_restricted,2,mean)

KBR$glycoprotein.score <- aaa[rownames(KBR@meta.data)]

#proteoglycan score
matrisome_set = read.table("https://raw.githubusercontent.com/mahmoudibrahim/KidneyMap/master/assets/public/ecm_genes_human.txt", sep = "\t", header = TRUE)
matrisome = matrisome_set$Gene.Symbol[which(matrisome_set$Category %in% c("Proteoglycans"))]
count.data = GetAssayData(KBR, slot = "data")
count.data = as.matrix(count.data[rownames(count.data) %in% toupper(matrisome) ,])
combined_matrix = rbind(count.data,average=apply(count.data,2,mean))
cor_matrix = cor(t(combined_matrix))
cor_vector = cor_matrix[,dim(cor_matrix)[1]]
count.data_restricted = count.data[rownames(count.data) %in% names(cor_vector[cor_vector >= 0.1]),]
aaa = apply(count.data_restricted,2,mean)

KBR$proteoglycan.score <- aaa[rownames(KBR@meta.data)]

KBR.sub <- subset(KBR, idents = 1:87)
VlnPlot(KBR.sub, features = c("aStr.score", "matrisome.score"), 
        ncol = 1, pt.size = 0, cols = cl.cols[levels(Idents(KBR.sub))]) + NoLegend()





###aEPI Marker Expression (clusters) (ED Fig 8e)

KBR.aEpi <- subset(KBR, idents = c(4:11,26:31,36:37))

Idents(KBR.aEpi) <- factor(Idents(KBR.aEpi), levels = c(4:11,31,36:37,26:27,30,28:29))

mark <- c(   "HAVCR1","VCAM1","PROM1","CST3","CD24",
             "CDHR5","MYO7B","USH1C", 
             "SLC7A8","SLC22A8","SLC22A6","SLC34A1","SLC22A5","SLC36A2","SLC5A12","SLC4A4","SLC13A3",
             "SLC12A6","SLC16A9",	"SLC6A13","SLC22A7","SLC5A10","SLC5A11","SLC6A12",
             "SLC28A1","SLC13A1","SLC16A12","SLC1A1","SLC17A1","SLC2A9","LRP2","CUBN",
             "KCNK13","WNK4",'CLCNKB',"SLC43A2","AKAP6","CPT1A","KCNJ3","SLC25A12",'MPC1',
             "KCNJ1","SLC9A3","SLC16A7","CLCNKA","KCNQ1","CLCN5","ACACB","SLC12A1","UMOD","EGF",
             "MET","LAMC1","ITGAV","NRP1","TNIK","PTPN1","SASH1","COL4A1","ACTN1","ITGA3",
             "DLGAP1","DLG1","CASK","FGFR1","FGFR1OP","FGFR2","FGF13","FGF14",  
             "TGFB1","TGFBR1","TGFBR2","SMAD1","SMAD2","SMAD3","SMAD4",
             "IL15","IL1R1","JAK1","STAT1","STAT3")

DotPlot(KBR.aEpi, features = mark,
        cols = c("light gray", "dodgerblue4"), dot.scale = 8) + RotatedAxis()





###snCv3/scCv3 Stromal Expression plot (Ext Data Fig 8f)
Idents(refquery) <- "subclass.l3"
RQ.Str <- subset(refquery, idents = c("FIB","dFIB","aFIB","MYOF","cycMYOF"))
Idents(RQ.Str) <- factor(Idents(RQ.Str), levels = c("FIB","dFIB","aFIB","MYOF","cycMYOF"))
ref.st.markers <- c("NEGR1","C7","FBLN5","LAMA2","MEG3","DCN",
                    "B2M","COL1A2","VIM","S100A6","AEBP1","COL1A1","SFRP1","CXCL12","PDGFRA",
                    "COL14A1",
                    "POSTN","FLRT2","FGF14","PRRX1","NAV3","IGF1","TWIST2",
                    "FAP","SULF1",
                    "ACTA2","COL6A3","COL16A1","GLI2","COL5A1","NTM","NOTCH3",
                    "MKI67","TOP2A")

DotPlot(RQ.Str, features = ref.st.markers,
        cols = c("light gray", "dodgerblue4"), dot.scale = 8) + RotatedAxis()





###SNARE2 RNA/AC - Degen State TFBS Activities (ED Fig 8b)
library(Signac)
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(swne)
library(chromfunks)
library(chromVAR)
library(viridis)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)

set.seed(1234)

##Prepare seurat objects
load("Kidney_Healthy-Injury_Cell_Atlas_SNARE2-RNA-AC_Seurat_03282022.rda")

#Combined snCv3/scCv3 Object
refquery <- merge(KBR, y = sc.KBR)
umap.coordinates <- rbind(Embeddings(object = KBR, reduction = "umap"),
                          Embeddings(object = sc.KBR, reduction = "ref.umap"))
refquery[["umap"]] <- CreateDimReducObject(embeddings = umap.coordinates, key = "umap_", assay = DefaultAssay(refquery))
refquery <- NormalizeData(refquery, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(refquery)
refquery <- ScaleData(refquery, features = all.genes, split.by = "id", assay = "RNA")
select.markers <- intersect(VariableFeatures(KBR), rownames(refquery))
VariableFeatures(refquery) <- select.markers


##Degenerative state TF activities
#PTs
Idents(KSAC) <- "subclass.l2"
Idents(KSAC) <- factor(Idents(KSAC), levels = c(
  "POD","dPOD","PEC","PT-S1","PT-S2","PT-S3","aPT","cycPT","dPT","dPT/DTL",
  "DTL2","DTL1","DTL3","dDTL3","ATL","dATL","aTAL1","aTAL2","M-TAL","dM-TAL",
  "C-TAL","dC-TAL","MD","DCT","dDCT","cycDCT","CNT","dCNT","cycCNT","C-PC",
  "M-PC","dM-PC","IMCD","dIMCD","C-IC-A","dC-IC-A","M-IC-A","tPC-IC","IC-B",
  "PapE","EC-GC","EC-AEA","EC-DVR","EC-PTC","dEC-PTC","EC-AVR","dEC","cycEC",
  "EC-LYM","MC","REN","VSMC","VSMC/P","dVSMC","MYOF","cycMYOF","FIB","M-FIB",
  "dM-FIB","aFIB","dFIB","B","PL","T","NKT","MAST","MAC-M2","cycMNP","MDC",
  "cDC","pDC","ncMON","N","SC/NEU"
))
KS.sub.pts13 <- subset(KSAC, idents = c("PT-S1","PT-S2","PT-S3",
                                        "dPT"))
KS.sub.pts13 <- RenameIdents(object = KS.sub.pts13, 
                             "PT-S1" = "PT", 
                             "PT-S2" = "PT",
                             "PT-S3" = "PT")
DefaultAssay(KS.sub.pts13) <- "chromvar"
tf.pts13.markers <- FindAllMarkers(
  object = KS.sub.pts13,
  only.pos = TRUE,
  test.use = 'LR',
  latent.vars = 'nCount_ATAC'
)


#TALs
KS.sub.tal <- subset(KSAC, idents = c("M-TAL","dM-TAL",
                                      "C-TAL","dC-TAL","MD"))
KS.sub.tal <- RenameIdents(object = KS.sub.tal, 
                           "M-TAL" = "TAL", 
                           "dM-TAL" = "dTAL",
                           "C-TAL" = "TAL", 
                           "dC-TAL" = "dTAL",
                           "MD" = "TAL")
DefaultAssay(KS.sub.tal) <- "chromvar"
tf.tal.markers <- FindAllMarkers(
  object = KS.sub.tal,
  only.pos = TRUE,
  test.use = 'LR',
  latent.vars = 'nCount_ATAC'
)


#CNT
KS.sub.cnt <- subset(KSAC, idents = c("CNT","dCNT"))
DefaultAssay(KS.sub.cnt) <- "chromvar"
tf.cnt.markers <- FindAllMarkers(
  object = KS.sub.cnt,
  only.pos = TRUE,
  test.use = 'LR',
  latent.vars = 'nCount_ATAC'
)

#PC
KS.sub.pc <- subset(KSAC, idents = c("C-PC","M-PC","dM-PC","IMCD","dIMCD"))
KS.sub.pc <- RenameIdents(object = KS.sub.pc, 
                          "C-PC" = "PC", 
                          "M-PC" = "PC",
                          "dM-PC" = "dPC", 
                          "IMCD" = "PC",
                          "dIMCD" = "dPC")
DefaultAssay(KS.sub.pc) <- "chromvar"
tf.pc.markers <- FindAllMarkers(
  object = KS.sub.pc,
  only.pos = TRUE,
  test.use = 'LR',
  latent.vars = 'nCount_ATAC'
)

#ICA
KS.sub.ic <- subset(KSAC, idents = c("C-IC-A",  "dC-IC-A",   "M-IC-A"))
KS.sub.ic <- RenameIdents(object = KS.sub.ic, 
                          "C-IC-A" = "IC", 
                          "M-IC-A" = "IC",
                          "dC-IC-A" = "dIC")
DefaultAssay(KS.sub.ic) <- "chromvar"
tf.ic.markers <- FindAllMarkers(
  object = KS.sub.ic,
  only.pos = TRUE,
  test.use = 'LR',
  latent.vars = 'nCount_ATAC'
)

#EC-PTC
KS.sub.ec <- subset(KSAC, idents = c("EC-PTC","dEC-PTC"))
DefaultAssay(KS.sub.ec) <- "chromvar"
tf.ec.markers <- FindAllMarkers(
  object = KS.sub.ec,
  only.pos = TRUE,
  test.use = 'LR',
  latent.vars = 'nCount_ATAC'
)

#FIB
KS.sub.fib <- subset(KSAC, idents = c("FIB","dFIB"))
DefaultAssay(KS.sub.fib) <- "chromvar"
tf.fib.markers <- FindAllMarkers(
  object = KS.sub.fib,
  only.pos = TRUE,
  test.use = 'LR',
  latent.vars = 'nCount_ATAC'
)

#Combine all tables
TF.list = list(tf.pts13.markers,
               tf.tal.markers,
               tf.cnt.markers,
               tf.pc.markers,
               tf.ic.markers,
               tf.ec.markers,
               tf.fib.markers)
celltype = 1:7

cons.TFs <- do.call(rbind, lapply(celltype, function(ct) {
  print(paste("Running for cell type:", ct))
  
  cl.mark <- TF.list[[ct]]
  cl.mark <- cl.mark[cl.mark$p_val < 0.05 & 
                       cl.mark$avg_log2FC > 0.35,]
  return(cl.mark)
  
}))


##Trim TFs to those expressed
Idents(refquery) <- "subclass.l1"
Idents(refquery) <- factor(Idents(refquery), levels = c(
  "POD","PEC","PT","DTL","ATL","TAL","DCT","CNT","PC","IC","PapE",
  "EC","VSM/P","FIB","IMM","NEU"))
refquery <- subset(refquery, idents = c("POD","PT","DTL","ATL","TAL","DCT","CNT","PC","IC",
                                        "EC","VSM/P","FIB"))
Idents(refquery) <- "state.l2"
refquery <- subset(refquery, idents = c("reference","degenerative"))
atac.rna.counts <- GetAssayData(object = refquery, slot = "counts")
atac.rna.idents <- Idents(object = refquery)

sn.cells <- rownames(refquery@meta.data[refquery@meta.data$id == "reference",])
sc.cells <- rownames(refquery@meta.data[refquery@meta.data$id == "query",])

sn.rna.avg.counts <- GetClusterAvg(atac.rna.counts[,sn.cells], atac.rna.idents[sn.cells], binarize = T)
sc.rna.avg.counts <- GetClusterAvg(atac.rna.counts[,sc.cells], atac.rna.idents[sc.cells], binarize = T)

hist(sn.rna.avg.counts)
quantile(sn.rna.avg.counts)
min.rna.cl.frac <- 0.2
sn.degen.expr <- sn.rna.avg.counts[,"degenerative"]
sn.degen.genes <- names(sn.degen.expr[sn.degen.expr > min.rna.cl.frac])
sc.degen.expr <- sc.rna.avg.counts[,"degenerative"]
sc.degen.genes <- names(sc.degen.expr[sc.degen.expr > min.rna.cl.frac])
degen.genes <- unique(c(sn.degen.genes,sc.degen.genes))

cons.TFs$motif <- cons.TFs$gene
cons.TFs$gene <- gsub("\\s*\\([^\\)]+\\)","",cons.TFs$gene)
cons.TFs$gene <- toupper(cons.TFs$gene)
cons.TFs$gene2 <- cons.TFs$gene

cons.TFs$gene2 <- gsub(":.*","",cons.TFs$gene)
cons.TFs$gene <- gsub(".*:","",cons.TFs$gene)

cons.TFs <- cons.TFs[cons.TFs$gene %in% degen.genes | cons.TFs$gene2 %in% degen.genes,]

inj.states <- c("dPT", "dTAL", "dCNT", "dPC", "dIC",
                "dEC-PTC", "dFIB")

inj.TF.counts <- data.frame(cons.TFs[cons.TFs$cluster %in% inj.states,] %>%
                              group_by(motif) %>%
                              tally())
inj.TF.counts <- inj.TF.counts[order(-inj.TF.counts$n),]

inf.TFs <- inj.TF.counts[inj.TF.counts$n >= 3, ]$motif

#Plot ref/degen states
Idents(KSAC) <- "subclass.l1"
Idents(KSAC) <- factor(Idents(KSAC), levels = c(
  "POD","PEC","PT","DTL","ATL","TAL","DCT","CNT","PC","IC","PapE",
  "EC","VSM/P","FIB","IMM","NEU"))
KSAC <- subset(KSAC, idents = c("POD","PT","DTL","ATL","TAL","DCT","CNT","PC","IC",
                                "EC","VSM/P","FIB"))
Idents(KSAC) <- "state.l2"
KSAC <- subset(KSAC, idents = c("reference","degenerative"))
Idents(KSAC) <- "subclass.l2"
Idents(KSAC) <- factor(Idents(KSAC), levels = c(
  "POD","dPOD","PEC","PT-S1","PT-S2","PT-S3","aPT","cycPT","dPT","dPT/DTL",
  "DTL2","DTL1","DTL3","dDTL3","ATL","dATL","aTAL1","aTAL2","M-TAL","dM-TAL",
  "C-TAL","dC-TAL","MD","DCT","dDCT","cycDCT","CNT","dCNT","cycCNT","C-PC",
  "M-PC","dM-PC","IMCD","dIMCD","C-IC-A","dC-IC-A","M-IC-A","tPC-IC","IC-B",
  "PapE","EC-GC","EC-AEA","EC-DVR","EC-PTC","dEC-PTC","EC-AVR","dEC","cycEC",
  "EC-LYM","MC","REN","VSMC","VSMC/P","dVSMC","MYOF","cycMYOF","FIB","M-FIB",
  "dM-FIB","aFIB","dFIB","B","PL","T","NKT","MAST","MAC-M2","cycMNP","MDC",
  "cDC","pDC","ncMON","N","SC/NEU"
))

DefaultAssay(KSAC) <- "chromvar"

DotPlot(KSAC, features = inf.TFs, dot.scale = 6, col.min = 0
) + RotatedAxis() + scale_color_viridis() + DarkTheme()





# SNARE2 - Cell type DARs and TFBS Activities (ED Fig 4) -----------------------

library(Signac)
library(Seurat)
library(dplyr)
library(Matrix)
library(methods)
library(swne)
library(ggplot2)
library(viridis)
library(chromfunks)
set.seed(1234)

###Prepare seurat objects
load("Kidney_Healthy-Injury_Cell_Atlas_SNARE2-RNA-AC_Seurat_03282022.rda")
load("color_factors.robj")



###Coverage Plots - ref-only clusters (ED Fig 4a)
DefaultAssay(KSAC) <- "ATAC"
Idents(KSAC) <- "subclass.l2"
Idents(KSAC) <- factor(Idents(KSAC), levels = c(
  "POD","dPOD","PEC","PT-S1","PT-S2","PT-S3","aPT","cycPT","dPT","dPT/DTL",
  "DTL2","DTL1","DTL3","dDTL3","ATL","dATL","aTAL1","aTAL2","M-TAL","dM-TAL",
  "C-TAL","dC-TAL","MD","DCT","dDCT","cycDCT","CNT","dCNT","cycCNT","C-PC",
  "M-PC","dM-PC","IMCD","dIMCD","C-IC-A","dC-IC-A","M-IC-A","tPC-IC","IC-B",
  "PapE","EC-GC","EC-AEA","EC-DVR","EC-PTC","dEC-PTC","EC-AVR","dEC","cycEC",
  "EC-LYM","MC","REN","VSMC","VSMC/P","dVSMC","MYOF","cycMYOF","FIB","M-FIB",
  "dM-FIB","aFIB","dFIB","B","PL","T","NKT","MAST","MAC-M2","cycMNP","MDC",
  "cDC","pDC","ncMON","N","SC/NEU"))
table(Idents(KSAC))
to.use <- c(
  "POD","PEC","PT-S1","PT-S2","PT-S3",
  "DTL2","DTL1","ATL",
  "M-TAL","C-TAL","MD","DCT",
  "CNT","C-PC","M-PC",
  "IMCD","C-IC-A","M-IC-A",
  "IC-B")



##Plot UMOD
umod.region <- LookupGeneCoords(KSAC, "UMOD")

cov_plot <- CoveragePlot(
  object = KSAC,
  idents = to.use,
  region = umod.region,
  annotation = FALSE,
  peaks = FALSE,
  links = FALSE)
cov_plot <- cov_plot+scale_fill_manual(values=as.character(sc.l2.cols[to.use]))

gene_plot <- AnnotationPlot(
  object = KSAC,
  region = umod.region
)

expr_plot <- ExpressionPlot(
  object = KSAC,
  idents = to.use,
  features = "UMOD",
  assay = "RNA"
)
expr_plot <- expr_plot+scale_fill_manual(values=as.character(sc.l2.cols[to.use]))

CombineTracks(
  plotlist = list(cov_plot, gene_plot),
  expression.plot = expr_plot,
  heights = c(10, 1),
  widths = c(10, 1)
)
#Repeat for SLC5A12, NPHS2, SH3GL3, SLC12A3, AQP2, SLC4A1






###subclass.l2 DEG-linked DARs (ED Fig 4b)
DefaultAssay(KSAC) <- "ATAC"
Idents(object = KSAC) <- "subclass.l2"

#Merge VSMC/P Populations
KSAC <- RenameIdents(KSAC, "VSMC" = "VSMC/P")
KSAC <- RenameIdents(KSAC, "MC" = "VSMC/P")
KSAC <- RenameIdents(KSAC, "REN" = "VSMC/P")

#Merge MD with C-TAL
KSAC <- RenameIdents(KSAC, "MD" = "C-TAL")

#Merge EC-AEA with EC-DVR
KSAC <- RenameIdents(KSAC, "EC-AEA" = "EC-DVR")


#Re-order Clusters
Idents(object = KSAC) <- factor(Idents(object = KSAC), levels = c(
  "POD","dPOD","PEC","PT-S1","PT-S2","PT-S3","aPT","cycPT","dPT","dPT/DTL",
  "DTL2","DTL1","DTL3","dDTL3","ATL","dATL","aTAL1","aTAL2","M-TAL","dM-TAL",
  "C-TAL","dC-TAL","MD","DCT","dDCT","cycDCT","CNT","dCNT","cycCNT","C-PC",
  "M-PC","dM-PC","IMCD","dIMCD","C-IC-A","dC-IC-A","M-IC-A","tPC-IC","IC-B",
  "PapE","EC-GC","EC-AEA","EC-DVR","EC-PTC","dEC-PTC","EC-AVR","dEC","cycEC",
  "EC-LYM","MC","REN","VSMC","VSMC/P","dVSMC","MYOF","cycMYOF","FIB","M-FIB",
  "dM-FIB","aFIB","dFIB","B","PL","T","NKT","MAST","MAC-M2","cycMNP","MDC",
  "cDC","pDC","ncMON","N","SC/NEU"
))

count.matrix <- GetAssayData(object = KSAC, slot = "counts")
count.matrix@x[count.matrix@x > 0] <- 1 #binarize counts

#Fix peak names
rownames(count.matrix) <- sub("-",":",rownames(count.matrix))

#Depth Corrected AC- Cluster DARs
##Call DARs against background of random cells with similar peak depth distribution
#Calculate the total reads in each cell based on accessibility matrix
#Estimate the distribution of total reads (depth distribution) for the cells belong to positive cluster
#Randomly sample cells from the rest of the clusters but in a weighted way so that we can sample the cells that have similar depth distribution in the positive cluster
#The weights for the background cells is estimated through the depth distribution so that we can know the probability of each cells in the distribution and the cells have similar depth have higher probability to be sampled
#Call DAR based on positive cells v.s. background cells

#' sample cells based on the distribution of depth in peak matrix
#' @param pmatrix peak matrix, rows are peak positions and columns are cells
#' @param cls.groups factor contains cell annotation for each celltype (cluster)
#' @param posCluster factor contains positive cluster
#' @param sampleSize number of sampled cells for negative class (background) - defaut equal to posCluster
#' @return vector of cells from positive and negative class
sampleCellDepth <- function(pmatrix, cls.groups, posCluster, sampleSize=length(posCluster)){
  depth <- Matrix::colSums(pmatrix)
  depthPos <- depth[names(depth) %in% names(posCluster)]
  depthNeg <- depth[-1*which(names(depth) %in% names(posCluster))]
  negCells <- names(cls.groups)[-1*which(names(cls.groups) %in% names(posCluster))]
  # fit distribution
  densityEst <- density(depthPos, kernel = "gaussian", bw = 1)
  weights <- approx(densityEst$x, densityEst$y, xout=depthNeg,
                    yright = 0.00001,
                    yleft = 0.00001)$y
  sampledCells <- negCells[sample(seq(length(negCells)), size = sampleSize,
                                  prob = weights, replace=FALSE)]
  totalCells <- c(names(posCluster), sampledCells)
  return(totalCells)
}
clusters.use <- levels(Idents(KSAC))
DAR.list <- lapply(clusters.use, function(cl) {
  print(paste("Running for Cluster:", cl))
  
  pmatrix = count.matrix
  cls.groups = Idents(object = KSAC)
  cls.x <- WhichCells(KSAC, idents = cl, downsample = 10000)
  posCluster = droplevels(cls.groups[cls.x])
  sampleSize = 10000
  cls.x.samp <- sampleCellDepth(pmatrix,cls.groups,posCluster,sampleSize)
  
  clusters <- cls.groups[cls.x.samp]
  current.cluster.ids <- levels(clusters)
  new.cluster.ids <- current.cluster.ids
  new.cluster.ids[!new.cluster.ids %in% levels(posCluster)] <- "other"
  clusters <- plyr::mapvalues(clusters, from = current.cluster.ids, to = new.cluster.ids)
  
  cls.x.DARs <- CalcDiffAccess(count.matrix[,cls.x.samp], clusters, min.p = 1e-100)
  
  cls.x.DARs[[cl]]
  
})
names(DAR.list) <- clusters.use
Top_acDARs <- lapply(DAR.list, function(df) subset(df, qval < 0.01))
sapply(Top_acDARs, nrow)

#Identify DEGs
KBR <- LoadH5Seurat("Kidney_Healthy-Injury_Cell_Atlas_snCv3_Seurat_03282022.h5Seurat")
Idents(object = KBR) <- "clusters"
Idents(object = KBR) <- factor(Idents(object = KBR), levels = 1:119)
KBR.sub <- subset(KBR, idents = 1:100)
Idents(KBR.sub) <- "subclass.l2"
Idents(KBR.sub) <- factor(Idents(KBR.sub), levels = c(
  "POD","dPOD","PEC","PT-S1","PT-S2","PT-S3","aPT","cycPT","dPT","dPT/DTL",
  "DTL2","DTL1","DTL3","dDTL3","ATL","dATL","aTAL1","aTAL2","M-TAL","dM-TAL",
  "C-TAL","dC-TAL","MD","DCT","dDCT","cycDCT","CNT","dCNT","cycCNT","C-PC",
  "M-PC","dM-PC","IMCD","dIMCD","C-IC-A","dC-IC-A","M-IC-A","tPC-IC","IC-B",
  "PapE","EC-GC","EC-AEA","EC-DVR","EC-PTC","dEC-PTC","EC-AVR","dEC","cycEC",
  "EC-LYM","MC","REN","VSMC","VSMC/P","dVSMC","MYOF","cycMYOF","FIB","M-FIB",
  "dM-FIB","aFIB","dFIB","B","PL","T","NKT","MAST","MAC-M2","cycMNP","MDC",
  "cDC","pDC","ncMON","N","SC/NEU"))
KBR.markers <- FindAllMarkers(KBR.sub, only.pos = TRUE, max.cells.per.ident = 500,
                              logfc.threshold = 0.25, min.pct = 0.25)

#Load Cicero connections (Calculated in "SNARE RNA/AC - Combined Seurat object")
conns.sub <- FilterConns(conns, min.coaccess = 0.1)
conns.sub$Peak1 <- sub("-",":",conns.sub$Peak1)
conns.sub$Peak2 <- sub("-",":",conns.sub$Peak2)

#Subset to DARs linked to marker genes
clusters.use <- names(which(sapply(Top_acDARs, nrow) > 100))

DAR.list.sub <- lapply(clusters.use, function(cl) {
  print(paste("Running for Cluster:", cl))
  
  peak.df <- Top_acDARs[[cl]]
  peaks <- rownames(peak.df)
  peaks <- peak2granges(peaks)
  
  genes <- KBR.markers[KBR.markers$cluster == cl,]$gene
  
  ## Link ATAC peaks to genes using cicero connections)
  peak.gene.df <- RegionGeneLinks(peaks, conns.sub, link.promoter = T,
                                  promoter.region = c(-3e3, 3e3),
                                  region.name = NULL, weight.col = "coaccess")
  
  peak.gene.df <- as.data.frame(peak.gene.df %>% group_by(region) %>% top_n(1, weight))
  peak.gene.df <- peak.gene.df[peak.gene.df$gene %in% genes, ]
  
  peak.df <- peak.df[peak.gene.df$region,]
  peak.df$gene <- peak.gene.df$gene
  peak.df$gene.weight <- peak.gene.df$weight
  peak.df
  
})
names(DAR.list.sub) <- clusters.use
Top_acDARs <- lapply(DAR.list.sub, function(df) subset(df, qval < 0.01 & logfc > 1))
sapply(Top_acDARs, nrow)
Top_acDARs.df <- do.call("rbind", lapply(Top_acDARs, as.data.frame)) 
cl<-unlist(lapply(rownames(Top_acDARs.df),function(x) unlist(strsplit(x,"[.]"))[1]))
loc <- unlist(lapply(rownames(Top_acDARs.df),function(x) unlist(strsplit(x,"[.]"))[2]))
Top_acDARs.df$cluster <- cl
Top_acDARs.df$site <- loc

#Visualization
cl.mark <- Top_acDARs.df
cl.mark$site <- gsub(":", "-", cl.mark$site)
cl.mark <- distinct(cl.mark, site, .keep_all = TRUE)
cl.mark %>% group_by(cluster) %>% top_n(10, logfc) -> top10

KSAC <- subset(KSAC, idents = unique(cl.mark$cluster))
ave.DAR <- AverageExpression(KSAC, assays = "ATAC", features = cl.mark$site, slot = "data" )

scaled <- t(scale(t(ave.DAR$ATAC)))
scaled <- scale(scaled)
scaled[scaled < 0] <- 0
scaled[scaled > 5] <- 5
ggHeat(t(scaled), rescaling = "none", clustering = "none", x.lab.size = 8, y.lab.size = 11,
) + scale_fill_viridis(option="magma")


#Motif enrichment for subclass.l2 DARs
DefaultAssay(KSAC) <- 'ATAC'
clusters <- names(which(table(cl.mark$cluster) > 2))
clusters <- levels(Idents(KSAC))[levels(Idents(KSAC)) %in% clusters]
enriched_motif_list <- lapply(clusters, function(cl) {
  EM <- FindMotifs(object = KSAC,
                   features = cl.mark[(cl.mark$cluster == cl),]$site)
  EM
})
names(enriched_motif_list) <- clusters

sc.em.df <- do.call("rbind", lapply(enriched_motif_list, as.data.frame)) 
cl<-unlist(lapply(rownames(sc.em.df),function(x) unlist(strsplit(x,"[.]"))[1]))
sc.em.df$cluster <- cl

MotifPlot(
  object = KSAC,
  motifs = c("Wt1", "HNF4A","THRB", "GRHL1", "ESRRB","GATA3","FOXI1","ETS1", "MEF2D","IRF1"),
  assay = 'ATAC'
)




###TFBS activities for expressed TFs (ED Fig 4c)
DefaultAssay(KSAC) <- 'chromvar'
tf.markers <- FindAllMarkers(
  object = KSAC,
  only.pos = TRUE,
  test.use = 'LR',
  latent.vars = 'nCount_ATAC'
)
DefaultAssay(KSAC) <- "RNA"
clusters <- levels(Idents(KSAC))
rna.counts <- GetAssayData(object = KSAC, slot = "counts")
rna.avg.counts <- GetClusterAvg(rna.counts, Idents(KSAC), binarize = T)
hist(rna.avg.counts)
quantile(rna.avg.counts)
min.rna.cl.frac <- 0.025
tf.markers$motif <- tf.markers$gene
tf.markers$gene <- gsub("\\s*\\([^\\)]+\\)","",tf.markers$gene) #remove " (var*)"
tf.markers$gene <- toupper(tf.markers$gene)

tf.markers.exp <- lapply(clusters, function(cl) {
  
  cl.expr <- rna.avg.counts[,cl]
  cl.expr <- names(cl.expr[cl.expr > min.rna.cl.frac])
  tf.markers.cl <- tf.markers[tf.markers$cluster == cl,]
  tf.markers.cl <- tf.markers.cl[tf.markers.cl$gene %in% cl.expr,]
  tf.markers.cl
})
names(tf.markers.exp) <- clusters

tf.markers.exp <- do.call("rbind", lapply(tf.markers.exp, as.data.frame)) 
rownames(tf.markers.exp) <- paste(tf.markers.exp$cluster,tf.markers.exp$motif, sep = ".")
tf.markers.exp <- na.omit(tf.markers.exp)

#Visualize
tf.mark <- tf.markers.exp[tf.markers.exp$p_val < 0.05,]
tf.mark <- distinct(tf.mark, gene, .keep_all = TRUE)
tf.mark %>% group_by(cluster) %>% top_n(5, avg_log2FC) -> top5 
DefaultAssay(KSAC) <- 'chromvar'
DotPlot(KSAC, features = top5$motif, dot.scale = 6, col.min = 0) + RotatedAxis() + scale_color_viridis(option="magma") + DarkTheme()




# snCv3/scCv3/SNARE - Expression/Regulation Dotplots (Fig 2h, Fig 5c, Fig 6c&g, ED Fig 11d, ED Fig 14e) ----------------------
library(Signac)
library(Seurat)
library(SeuratDisk)
library(dplyr)
library(Matrix)
library(ggplot2)
library(swne)
library(viridis)
library(chromfunks)

set.seed(1234)
load("color_factors.robj")
load("sc-sn_int_color_factors.robj")

###Prepare seurat objects
load("Kidney_Healthy-Injury_Cell_Atlas_SNARE2-RNA-AC_Seurat_03282022.rda")

KBR <- LoadH5Seurat("Kidney_Healthy-Injury_Cell_Atlas_snCv3_Seurat_03282022.h5Seurat")
Idents(object = KBR) <- "clusters"
Idents(object = KBR) <- factor(Idents(object = KBR), levels = 1:119)
KBR <- subset(KBR, idents = 1:100)
KBR$id <- "ref"
sc.KBR <- LoadH5Seurat("Kidney_Healthy-Injury_Cell_Atlas_scCv3_Seurat_03282022.h5Seurat")
colnames(sc.KBR@meta.data)[colnames(sc.KBR@meta.data) == "state"] <- "state.l1"

#Combined snCv3/scCv3 Object
refquery <- merge(KBR, y = sc.KBR)
umap.coordinates <- rbind(Embeddings(object = KBR, reduction = "umap"),
                          Embeddings(object = sc.KBR, reduction = "ref.umap"))
refquery[["umap"]] <- CreateDimReducObject(embeddings = umap.coordinates, key = "umap_", assay = DefaultAssay(refquery))
refquery <- NormalizeData(refquery, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(refquery)
refquery <- ScaleData(refquery, features = all.genes, split.by = "id", assay = "RNA")
select.markers <- intersect(VariableFeatures(KBR), rownames(refquery))
VariableFeatures(refquery) <- select.markers
Idents(refquery) <- "subclass.l3"
#Merge to concensus labels 
refquery <- RenameIdents(refquery, 
                         'PT-S1' = 'PT-S1/2',
                         'PT-S2' = 'PT-S1/2',
                         'NKT' = 'NKC/T',
                         'dPT/DTL' = 'dPT')
refquery$subclass.l3 <- Idents(refquery)



###Renal corpuscle dotplot (Fig. 2h)
RQ.GLOM <- subset(refquery, idents = c(
  "POD","PEC","PT-S1/2","MC","EC-GC","EC-AEA","MD","REN","VSMC"))
Idents(RQ.GLOM) <- factor(Idents(RQ.GLOM), levels = c(
  "POD","EC-GC","MC","PEC","PT-S1/2","EC-AEA","VSMC","REN","MD"))

DotPlot(RQ.GLOM, features = c("NPHS2","EMCN","POSTN","VCAM1","SLC5A12","PALMD","MYH11","REN","NOS1"),
        cols = c("light gray", "dodgerblue4"), dot.scale = 8) + RotatedAxis()
rm(RQ.GLOM)


###aSTR dotplot (ED Fig 8f)
RQ.Str <- subset(refquery, idents = c("FIB","dFIB","aFIB","MYOF","cycMYOF"))
Idents(RQ.Str) <- factor(Idents(RQ.Str), levels = c("FIB","dFIB","aFIB","MYOF","cycMYOF"))
ref.st.markers <- c("NEGR1","C7","FBLN5","LAMA2","MEG3","DCN",                                 
                    "B2M","COL1A2","VIM","S100A6","AEBP1","COL1A1","SFRP1","CXCL12","PDGFRA",               #FIB
                    "COL14A1",
                    "POSTN","FLRT2","FGF14","PRRX1","NAV3","IGF1","TWIST2",                    
                    "FAP","SULF1",                                                                           
                    "ACTA2","COL6A3","COL16A1","GLI2","COL5A1","NTM","NOTCH3",
                    "MKI67","TOP2A")

DotPlot(RQ.Str, features = ref.st.markers,
        cols = c("light gray", "dodgerblue4"), dot.scale = 8) + RotatedAxis()
rm(RQ.Str)





###Aging/SASP Gene set scores (Fig. 6g) 
#Prepare objects for score calculation
ref <- KBR
query <- sc.KBR
order <- c(
  "POD","dPOD","PEC","PT-S1/2","PT-S3","aPT","cycPT","dPT",
  "dPT/DTL","DTL2","DTL1","DTL3","dDTL3","ATL","dATL","aTAL1",
  "aTAL2","M-TAL","dM-TAL","C-TAL","dC-TAL","MD","DCT1","DCT2",
  "dDCT","cycDCT","CNT","CNT-PC","dCNT","cycCNT","CCD-PC","OMCD-PC",
  "dOMCD-PC","IMCD","dIMCD","CCD-IC-A","CNT-IC-A","dC-IC-A","OMCD-IC-A",
  "tPC-IC","IC-B","PapE","EC-GC","EC-AEA","EC-DVR","EC-PTC","dEC-PTC",
  "EC-AVR","dEC","cycEC","EC-LYM","MC","REN","VSMC","VSMC/P","dVSMC",
  "MYOF","cycMYOF","FIB","M-FIB","dM-FIB","aFIB","dFIB","B","PL","T",
  "NKC/T","cycNKC/T","MAST","MAC-M2","cycMNP","MDC","cDC","pDC","ncMON","N","SC/NEU")

sen.genes <- read.delim("Scenescence_Marker_Gene_Sets.txt",sep="\t",header=TRUE)

##Tabula.Muris
count.data = GetAssayData(ref, slot = "data")
sen = rownames(count.data)[rownames(count.data) %in% toupper(sen.genes$Tabula.Muris)]
count.data = as.matrix(count.data[rownames(count.data) %in% sen ,])
combined_matrix = rbind(count.data,average=apply(count.data,2,mean))
cor_matrix = cor(t(combined_matrix))
cor_vector = cor_matrix[,dim(cor_matrix)[1]]
count.data_restricted = count.data[rownames(count.data) %in% names(cor_vector[cor_vector >= 0.1]),]
aaa = apply(count.data_restricted,2,mean)
ref$tabula.muris.score <- aaa[rownames(ref@meta.data)]

count.data = GetAssayData(query, slot = "data")
sen = rownames(count.data)[rownames(count.data) %in% toupper(sen.genes$Tabula.Muris)]
count.data = as.matrix(count.data[rownames(count.data) %in% sen ,])
combined_matrix = rbind(count.data,average=apply(count.data,2,mean))
cor_matrix = cor(t(combined_matrix))
cor_vector = cor_matrix[,dim(cor_matrix)[1]]
count.data_restricted = count.data[rownames(count.data) %in% names(cor_vector[cor_vector >= 0.1]),]
aaa = apply(count.data_restricted,2,mean)
query$tabula.muris.score <- aaa[rownames(query@meta.data)]

##Takemon
count.data = GetAssayData(ref, slot = "data")
sen = rownames(count.data)[rownames(count.data) %in% toupper(sen.genes$Takemon)]
count.data = as.matrix(count.data[rownames(count.data) %in% sen ,])
combined_matrix = rbind(count.data,average=apply(count.data,2,mean))
cor_matrix = cor(t(combined_matrix))
cor_vector = cor_matrix[,dim(cor_matrix)[1]]
count.data_restricted = count.data[rownames(count.data) %in% names(cor_vector[cor_vector >= 0.1]),]
aaa = apply(count.data_restricted,2,mean)
ref$Takemon.score <- aaa[rownames(ref@meta.data)]

count.data = GetAssayData(query, slot = "data")
sen = rownames(count.data)[rownames(count.data) %in% toupper(sen.genes$Takemon)]
count.data = as.matrix(count.data[rownames(count.data) %in% sen ,])
combined_matrix = rbind(count.data,average=apply(count.data,2,mean))
cor_matrix = cor(t(combined_matrix))
cor_vector = cor_matrix[,dim(cor_matrix)[1]]
count.data_restricted = count.data[rownames(count.data) %in% names(cor_vector[cor_vector >= 0.1]),]
aaa = apply(count.data_restricted,2,mean)
query$Takemon.score <- aaa[rownames(query@meta.data)]

##Ruscetti
count.data = GetAssayData(ref, slot = "data")
sen = rownames(count.data)[rownames(count.data) %in% toupper(sen.genes$Ruscetti)]
count.data = as.matrix(count.data[rownames(count.data) %in% sen ,])
combined_matrix = rbind(count.data,average=apply(count.data,2,mean))
cor_matrix = cor(t(combined_matrix))
cor_vector = cor_matrix[,dim(cor_matrix)[1]]
count.data_restricted = count.data[rownames(count.data) %in% names(cor_vector[cor_vector >= 0.1]),]
aaa = apply(count.data_restricted,2,mean)
ref$Ruscetti.score <- aaa[rownames(ref@meta.data)]

count.data = GetAssayData(query, slot = "data")
sen = rownames(count.data)[rownames(count.data) %in% toupper(sen.genes$Ruscetti)]
count.data = as.matrix(count.data[rownames(count.data) %in% sen ,])
combined_matrix = rbind(count.data,average=apply(count.data,2,mean))
cor_matrix = cor(t(combined_matrix))
cor_vector = cor_matrix[,dim(cor_matrix)[1]]
count.data_restricted = count.data[rownames(count.data) %in% names(cor_vector[cor_vector >= 0.1]),]
aaa = apply(count.data_restricted,2,mean)
query$Ruscetti.score <- aaa[rownames(query@meta.data)]

##Basisty
count.data = GetAssayData(ref, slot = "data")
sen = rownames(count.data)[rownames(count.data) %in% toupper(sen.genes$Basisty)]
count.data = as.matrix(count.data[rownames(count.data) %in% sen ,])
combined_matrix = rbind(count.data,average=apply(count.data,2,mean))
cor_matrix = cor(t(combined_matrix))
cor_vector = cor_matrix[,dim(cor_matrix)[1]]
count.data_restricted = count.data[rownames(count.data) %in% names(cor_vector[cor_vector >= 0.1]),]
aaa = apply(count.data_restricted,2,mean)
ref$Basisty.score <- aaa[rownames(ref@meta.data)]

count.data = GetAssayData(query, slot = "data")
sen = rownames(count.data)[rownames(count.data) %in% toupper(sen.genes$Basisty)]
count.data = as.matrix(count.data[rownames(count.data) %in% sen ,])
combined_matrix = rbind(count.data,average=apply(count.data,2,mean))
cor_matrix = cor(t(combined_matrix))
cor_vector = cor_matrix[,dim(cor_matrix)[1]]
count.data_restricted = count.data[rownames(count.data) %in% names(cor_vector[cor_vector >= 0.1]),]
aaa = apply(count.data_restricted,2,mean)
query$Basisty.score <- aaa[rownames(query@meta.data)]

##Add to main Refquery object
refquery@meta.data <- refquery@meta.data[c(rownames(ref@meta.data),rownames(query@meta.data)),]
refquery@meta.data$tabula.muris.score <- c(ref@meta.data$tabula.muris.score,query@meta.data$tabula.muris.score)
refquery@meta.data$Takemon.score <- c(ref@meta.data$Takemon.score,query@meta.data$Takemon.score)
refquery@meta.data$Ruscetti.score <- c(ref@meta.data$Ruscetti.score,query@meta.data$Ruscetti.score)
refquery@meta.data$Basisty.score <- c(ref@meta.data$Basisty.score,query@meta.data$Basisty.score)

##Score Plots (ED Fig 14e)
Idents(refquery) <- "subclass.l3"
Idents(refquery) <- factor(Idents(refquery), levels = order)
VlnPlot(refquery, features = c("tabula.muris.score","Takemon.score","Ruscetti.score","Basisty.score"), 
        ncol = 1, pt.size = 0, cols = hsc.l3.cols[levels(Idents(refquery))]) + NoLegend()






###TFBS activities for aPT states (Supplementary Table 20, ED Fig 11d)
mod.assign <- readRDS("pt.snare.module.assignment.rds")
KSAC@meta.data$module <- mod.assign[rownames(KSAC@meta.data)]

Idents(KSAC) <- "module"
KSAC.PT.as <- subset(KSAC, idents = c("black","red","magenta","yellow","turquoise"))
Idents(KSAC.PT.as) <- factor(Idents(KSAC.PT.as), levels = c("black","red","magenta","yellow","turquoise"))

DefaultAssay(KSAC.PT.as) <- "chromvar"
tf.markers <- FindAllMarkers(
  object = KSAC.PT.as,
  max.cells.per.ident = 500,
  only.pos = TRUE,
  test.use = 'LR',
  latent.vars = 'nCount_ATAC'
)

#Subset to TFs enriched and expressed (snCv3/scCv3) 
sn.pt.mod <- readRDS("pt.sn.module.assignment.rds")
sc.pt.mod <- readRDS("pt.sc.module.assignment.rds")
refquery@meta.data[names(sn.pt.mod),]$pt.module <- sn.pt.mod
refquery@meta.data[names(sc.pt.mod),]$pt.module <- sc.pt.mod
Idents(refquery) <- "pt.module"
RQ.PT <- subset(refquery, idents = c("black","red","magenta","yellow","turquoise"))
Idents(RQ.PT) <- factor(Idents(RQ.PT), levels = c("black","red","magenta","yellow","turquoise"))

DefaultAssay(RQ.PT) <- "RNA"
clusters <- levels(Idents(RQ.PT))
rna.counts <- GetAssayData(object = RQ.PT, slot = "counts")
rna.avg.counts <- GetClusterAvg(rna.counts, Idents(RQ.PT), binarize = T)
hist(rna.avg.counts)
quantile(rna.avg.counts)
min.rna.cl.frac <- 0.025
tf.markers$motif <- tf.markers$gene
tf.markers$gene <- gsub("\\s*\\([^\\)]+\\)","",tf.markers$gene) #remove " (var*)"
tf.markers$gene <- toupper(tf.markers$gene)
tf.markers$gene2 <- tf.markers$gene

tf.markers$gene2 <- gsub(":.*","",tf.markers$gene)
tf.markers$gene <- gsub(".*:","",tf.markers$gene)

tf.markers.exp <- lapply(clusters, function(cl) {
  
  cl.expr <- rna.avg.counts[,cl]
  cl.expr <- names(cl.expr[cl.expr > min.rna.cl.frac])
  tf.markers.cl <- tf.markers[tf.markers$cluster == cl,]
  tf.markers.cl <- tf.markers.cl[tf.markers.cl$gene %in% cl.expr | tf.markers.cl$gene2 %in% cl.expr,]
  tf.markers.cl
})
names(tf.markers.exp) <- clusters
tf.markers.exp <- do.call("rbind", lapply(tf.markers.exp, as.data.frame)) 
rownames(tf.markers.exp) <- paste(tf.markers.exp$cluster,tf.markers.exp$motif, sep = ".")
tf.markers.exp

#Plot representative gene expression and TF activities
gene.sub <- c("HNF4A", "THRB", "NR3C2", "MAF", "NFIC", 
              "JUND","JUN","FOSL2","JUNB","FOS","BACH2","ATF2",
              "TGFBR1","TGFBR2","SMAD3","SMAD2",
              "NFKB1","RELB","REL","ATF3","NFKB2","RELA",
              "JAK1","STAT1","STAT2",
              "CREB5","CREB1","CREM",
              "ELF3","ETV6","PAX2","ETS1","EGR1","MEIS1")

DotPlot(RQ.PT, features = gene.sub,
        cols = c("light gray", "dodgerblue4"), dot.scale = 8) + RotatedAxis()

TFBS.sub <- c("HNF4A(var.2)","THRB","NR3C2","MAF","NFIC",
              "JUND(var.2)","JUN","FOSL2","FOS::JUNB","FOS","BACH2","ATF2",
              "Smad2::Smad3",
              "NFKB1","RELB","REL","ATF3","NFKB2","RELA",
              "STAT1::STAT2",
              "Creb5","CREB1","CREM",
              "ELF3","ETV6","Pax2","ETS1","EGR1","MEIS1(var.2)"
)

DefaultAssay(KSAC.PT.as) <- "chromvar"
DotPlot(KSAC.PT.as, features = TFBS.sub, dot.scale = 6, col.min = 0
) + RotatedAxis() + scale_color_viridis() + DarkTheme()







###TFBS activities for aTAL states (Supplementary Table 20, Fig 5c)
snare.mod.assign <- readRDS("tal.snare.module.assignment.rds")
snare.mod.assign <- snare.mod.assign[rownames(KSAC@meta.data)]
KSAC@meta.data$tal.module <- snare.mod.assign

Idents(KSAC) <- "tal.module"
KSAC.TAL <- subset(KSAC, idents = c("black", "pink", "brown", "yellow", "blue"))
Idents(KSAC.TAL) <- factor(Idents(KSAC.TAL), levels = c("black", "pink", "brown", "yellow", "blue"))

DefaultAssay(KSAC.TAL) <- "chromvar"
tf.markers <- FindAllMarkers(
  object = KSAC.TAL,
  max.cells.per.ident = 500,
  only.pos = TRUE,
  test.use = 'LR',
  latent.vars = 'nCount_ATAC'
)

#Subset to TFs enriched and expressed (snCv3/scCv3) 
sn.sc.mod.assign <- readRDS("tal.sn.sc.module.assignment.rds")
sn.sc.mod.assign <- sn.sc.mod.assign[rownames(refquery@meta.data)]

refquery@meta.data$tal.module <- sn.sc.mod.assign
Idents(refquery) <- "tal.module"
RQ.TAL <- subset(refquery, idents = c("black", "pink", "brown", "yellow", "blue"))
Idents(RQ.TAL) <- factor(Idents(RQ.TAL), levels = c("black", "pink", "brown", "yellow", "blue"))

DefaultAssay(RQ.TAL) <- "RNA"
clusters <- levels(Idents(RQ.TAL))
rna.counts <- GetAssayData(object = RQ.TAL, slot = "counts")
rna.avg.counts <- GetClusterAvg(rna.counts, Idents(RQ.TAL), binarize = T)
hist(rna.avg.counts)
quantile(rna.avg.counts)
min.rna.cl.frac <- 0.025
tf.markers$motif <- tf.markers$gene
tf.markers$gene <- gsub("\\s*\\([^\\)]+\\)","",tf.markers$gene) #remove " (var*)"
tf.markers$gene <- toupper(tf.markers$gene)
tf.markers$gene2 <- tf.markers$gene

tf.markers$gene2 <- gsub(":.*","",tf.markers$gene)
tf.markers$gene <- gsub(".*:","",tf.markers$gene)

tf.markers.exp <- lapply(clusters, function(cl) {
  
  cl.expr <- rna.avg.counts[,cl]
  cl.expr <- names(cl.expr[cl.expr > min.rna.cl.frac])
  tf.markers.cl <- tf.markers[tf.markers$cluster == cl,]
  tf.markers.cl <- tf.markers.cl[tf.markers.cl$gene %in% cl.expr | tf.markers.cl$gene2 %in% cl.expr,]
  tf.markers.cl
})
names(tf.markers.exp) <- clusters
tf.markers.exp <- do.call("rbind", lapply(tf.markers.exp, as.data.frame)) 
rownames(tf.markers.exp) <- paste(tf.markers.exp$cluster,tf.markers.exp$motif, sep = ".")
tf.markers.exp


#Plot Representative Subset - RNA and TFBS activities
gene.sub <- c("ESRRG","PPARGC1A","ESRRB","ESRRA",
              "JUND","FOS","JUN","JUNB","BACH1","FOSL2","BACH2","ATF2",
              "ARHGAP29","ROCK2","RAC1","CDC42","RHOA","RHOB","RHOC",
              "STAT3","JAK1","STAT5B","STAT1",
              "SMAD3","SMAD2","TGFBR2","TGFBR1","ITGB6","ITGB8","ITGAV","SKIL",
              "NFKBIZ","NFKB1","RELB","REL",
              "BICC1","DCDC2","TCF7L2","TNIK","PRICKLE1",
              "NR2F2","NR3C1","TFAP2B","KCTD1","ETV6","ZEB1","EHF","ELF3","ETS1","RFX2","CEBPD",
              "FGF13","FGFR1","FGF2",
              "EGF","ERBB4","EGFR","NRG1","ERBB2","PLSCR1","NRG3",
              "MET","MACC1"
              
)
DotPlot(RQ.TAL, features = gene.sub,
        cols = c("light gray", "dodgerblue4"), dot.scale = 8) + RotatedAxis()

TFBS.sub <- c("Esrrg","ESRRB","ESRRA",
              "FOS::JUND","FOSL1::JUN","FOSB::JUNB","BACH1","FOSL2","BACH2","ATF2",
              
              "STAT3","Stat5a::Stat5b","STAT1",
              "Smad2::Smad3",
              
              "REL","RELB","NFKB1",
              "TCF7L2",
              "NR2F2","NR3C1","TFAP2B","ETV6","ZEB1","EHF","ELF3","ETS1","RFX2","CEBPD"
)

DefaultAssay(KSAC.TAL) <- "chromvar"
DotPlot(KSAC.TAL, features = TFBS.sub, dot.scale = 6, col.min = 0
) + RotatedAxis() + scale_color_viridis() + DarkTheme()




###Ligand Receptor Analysis: TAL genes (Fig. 6c)
Idents(refquery) <- "tal.module"
RQ.TAL <- subset(refquery, idents = c("pink", "brown", "yellow", "blue"))
Idents(RQ.TAL) <- factor(Idents(RQ.TAL), levels = c("pink", "brown", "yellow", "blue"))

DotPlot(RQ.TAL, features = c("SLC12A1","EGF","IGF1R","NECTIN2","NRG1","NRG3","BMP6","NRP1","SEMA3C"),
        cols = c("light gray", "dodgerblue4"), dot.scale = 8) + RotatedAxis()




###Ligand Receptor Analysis: STR genes (Fig. 6c)
Idents(refquery) <- "subclass.l3"
RQ.STR <- subset(refquery, idents = c("FIB","aFIB","MYOF"))
Idents(RQ.STR) <- factor(Idents(RQ.STR), levels = c("FIB","aFIB","MYOF"))
DotPlot(RQ.STR, features = c("FBLN5","FLRT2","FGF14","IGF1","IGFBP6","IGFBP5","ERBB4","BMPR2","BMPR1A","ACVR1","BMPR1B",
                             "NRP1","SEMA3C"), 
        cols = c("light gray", "dodgerblue4"), dot.scale = 8) + RotatedAxis()




###Ligand Receptor Analysis: IMM genes (Fig. 6c)
Idents(refquery) <- "subclass.l3"
RQ.IMM <- subset(refquery, idents = c("PL","T","NKC/T","MAC-M2","MDC","ncMON"))
Idents(RQ.IMM) <- factor(Idents(RQ.IMM), levels = c("PL","T","NKC/T","MAC-M2","MDC","ncMON"))
DotPlot(RQ.IMM, features = c("CD226","ERBB4","NRP1"), 
        cols = c("light gray", "dodgerblue4"), dot.scale = 8) + RotatedAxis()



###Ligand Receptor Analysis: TAL genes (Fig. 6c)
RQ.TAL <- subset(RQ.TAL, idents = c("pink", "brown", "yellow", "blue"))
Idents(RQ.TAL) <- factor(Idents(RQ.TAL), levels = c("pink", "brown", "yellow", "blue"))

DotPlot(RQ.TAL, features = c("SLC12A1","EGF","IGF1R","NECTIN2","NRG1","NRG3","BMP6","NRP1","SEMA3C"),
        cols = c("light gray", "dodgerblue4"), dot.scale = 8) + RotatedAxis()




###Aging and SASP expresssion (Fig. 6g)
#PT Modules
sn.pt.mod <- readRDS("pt.sn.module.assignment.rds")
sc.pt.mod <- readRDS("pt.sc.module.assignment.rds")
pt.mod <- c(sn.pt.mod,sc.pt.mod)
cells <- names(pt.mod)
pt.mod <- paste0("PT_", pt.mod)
names(pt.mod) <- cells
refquery@meta.data$module <- NA
refquery@meta.data[names(pt.mod),]$module <- pt.mod

#TAL modules
tal.mod <- readRDS("tal.sn.sc.module.assignment.rds")
cells <- names(tal.mod)
tal.mod <- paste0("TAL_", tal.mod)
names(tal.mod) <- cells
refquery@meta.data[names(tal.mod),]$module <- tal.mod

Idents(refquery) <- "module"
RQ.aEpi <- subset(refquery, idents = c("PT_black","PT_red","PT_magenta","PT_yellow","PT_turquoise",
                                       "TAL_black","TAL_pink","TAL_brown","TAL_yellow","TAL_blue"))
Idents(RQ.aEpi) <- factor(Idents(RQ.aEpi), levels = c("PT_black","PT_red","PT_magenta","PT_yellow","PT_turquoise",
                                                      "TAL_black","TAL_pink","TAL_brown","TAL_yellow","TAL_blue"))
senescent.genes <- c("CDKN1A", "CDKN1B", "CDKN2A")
aEpi <- c("CDH6","HAVCR1","VCAM1","PROM1","ITGB6","CCL2")

DotPlot(RQ.aEpi, features = aEpi,
        cols = c("light gray", "dodgerblue4"), dot.scale = 8) + RotatedAxis()
DotPlot(RQ.aEpi, features = senescent.genes,
        cols = c("light gray", "#99000d"), dot.scale = 8) + RotatedAxis()

VlnPlot(RQ.aEpi, features = c("tabula.muris.score","Takemon.score","Ruscetti.score","Basisty.score"), 
        ncol = 1, pt.size = 0, cols = c("black","red","magenta","yellow","turquoise",
                                        "black","pink","brown","yellow","blue")) + NoLegend()






###Stress TF Activities (Fig. 6g)
#PT Modules
pt.mod <- readRDS("pt.snare.module.assignment.rds")
cells <- names(pt.mod)
pt.mod <- paste0("PT_", pt.mod)
names(pt.mod) <- cells
KSAC@meta.data$module <- NA
KSAC@meta.data[names(pt.mod),]$module <- pt.mod

#TAL modules
tal.mod <- readRDS("tal.snare.module.assignment.rds")
cells <- names(tal.mod)
tal.mod <- paste0("TAL_", tal.mod)
names(tal.mod) <- cells
KSAC@meta.data[names(tal.mod),]$module <- tal.mod

Idents(KSAC) <- "module"
KSAC.as <- subset(KSAC, idents = c("PT_black","PT_red","PT_magenta","PT_yellow","PT_turquoise",
                                   "TAL_black", "TAL_pink", "TAL_brown", "TAL_yellow", "TAL_blue"))
Idents(KSAC.as) <- factor(Idents(KSAC.as), levels = c("PT_black","PT_red","PT_magenta","PT_yellow","PT_turquoise",
                                                      "TAL_black", "TAL_pink", "TAL_brown", "TAL_yellow", "TAL_blue"))
TFBS.sub <- c("XBP1","ATF4","ATF6","HIF1A")

DefaultAssay(KSAC.as) <- "chromvar"
DotPlot(KSAC.as, features = TFBS.sub, dot.scale = 6, col.min = 0
) + RotatedAxis() + scale_color_viridis() + DarkTheme()

##



# SNARE2 - GWAS and eQTL Analyses (Fig. 6f, ED Fig 14c-d) ---------------------------------------------------------

library(chromVAR)
library(gchromVAR)
library(SummarizedExperiment)
library(data.table)
library(BiocParallel)
library(BSgenome.Hsapiens.UCSC.hg38)
library(chromfunks)
library(swne)
library(ggplot2)
library(Seurat)
library(Signac)
library(Matrix)
library(dplyr)
library(cowplot)
library(ChIPseeker)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(viridis)
library(chromVARmotifs)
library(motifmatchr)
library(igraph)

## Global parameters
register(MulticoreParam(6))
n.runs <- 20


###Prepare seurat objects
##Combined snCv3/scCv3 Object
KBR <- LoadH5Seurat("Kidney_Healthy-Injury_Cell_Atlas_snCv3_Seurat_03282022.h5Seurat")
Idents(object = KBR) <- "clusters"
Idents(object = KBR) <- factor(Idents(object = KBR), levels = 1:119)
KBR <- subset(KBR, idents = 1:100)
KBR$id <- "ref"
sc.KBR <- LoadH5Seurat("Kidney_Healthy-Injury_Cell_Atlas_scCv3_Seurat_03282022.h5Seurat")
colnames(sc.KBR@meta.data)[colnames(sc.KBR@meta.data) == "state"] <- "state.l1"
refquery <- merge(KBR, y = sc.KBR)
umap.coordinates <- rbind(Embeddings(object = KBR, reduction = "umap"),
                          Embeddings(object = sc.KBR, reduction = "ref.umap"))
refquery[["umap"]] <- CreateDimReducObject(embeddings = umap.coordinates, key = "umap_", assay = DefaultAssay(refquery))
refquery <- NormalizeData(refquery, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(refquery)
refquery <- ScaleData(refquery, features = all.genes, split.by = "id", assay = "RNA")
select.markers <- intersect(VariableFeatures(KBR), rownames(refquery))
VariableFeatures(refquery) <- select.markers
Idents(refquery) <- "subclass.l3"
#Merge to concensus labels 
refquery <- RenameIdents(refquery, 
                         'PT-S1' = 'PT-S1/2',
                         'PT-S2' = 'PT-S1/2',
                         'NKT' = 'NKC/T',
                         'dPT/DTL' = 'dPT')
refquery$subclass.l3 <- Idents(refquery)

##SNARE Object
load("Kidney_Healthy-Injury_Cell_Atlas_SNARE2-RNA-AC_Seurat_03282022.rda")
DefaultAssay(KSAC) <- "ATAC"
Idents(object = KSAC) <- "subclass.l2"

#Merge VSMC/P Populations
KSAC <- RenameIdents(KSAC, "VSMC" = "VSMC/P")
KSAC <- RenameIdents(KSAC, "MC" = "VSMC/P")
KSAC <- RenameIdents(KSAC, "REN" = "VSMC/P")

#Merge MD with C-TAL
KSAC <- RenameIdents(KSAC, "MD" = "C-TAL")

#Merge EC-AEA with EC-DVR
KSAC <- RenameIdents(KSAC, "EC-AEA" = "EC-DVR")

#Re-order Clusters
Idents(object = KSAC) <- factor(Idents(object = KSAC), levels = c(
  "POD","dPOD","PEC","PT-S1","PT-S2","PT-S3","aPT","cycPT","dPT","dPT/DTL",
  "DTL2","DTL1","DTL3","dDTL3","ATL","dATL","aTAL1","aTAL2","M-TAL","dM-TAL",
  "C-TAL","dC-TAL","MD","DCT","dDCT","cycDCT","CNT","dCNT","cycCNT","C-PC",
  "M-PC","dM-PC","IMCD","dIMCD","C-IC-A","dC-IC-A","M-IC-A","tPC-IC","IC-B",
  "PapE","EC-GC","EC-AEA","EC-DVR","EC-PTC","dEC-PTC","EC-AVR","dEC","cycEC",
  "EC-LYM","MC","REN","VSMC","VSMC/P","dVSMC","MYOF","cycMYOF","FIB","M-FIB",
  "dM-FIB","aFIB","dFIB","B","PL","T","NKT","MAST","MAC-M2","cycMNP","MDC",
  "cDC","pDC","ncMON","N","SC/NEU"
))

#Subset to cell types with at least 100 DARs
#DAR.list from "SNARE2 - Cell type DARs and TFBS Activities"
Top_acDARs <- lapply(DAR.list, function(df) subset(df, qval < 0.01))
sapply(Top_acDARs, nrow)
clusters.use <- names(which(sapply(Top_acDARs, nrow) > 100))
KSAC <- subset(KSAC, idents = clusters.use)

##Assess causal variant enrichments
sc.counts <- GetAssayData(KSAC, assay = "ATAC")
sc.counts@x[sc.counts@x > 1] <- 1 ## Binarize the data
dim(sc.counts)

clusters <- Idents(KSAC)
clusters <- droplevels(clusters[!is.na(clusters)])

counts <- getPseudobulk(sc.counts[,names(clusters)], clusters)
counts <- counts[rowSums(counts) > 0,]
peaks <- peak2granges(rownames(counts), delim = c("-", "-"))

## Clean up unneeded objects
rm(sc.counts); invisible(gc())

## Build summarized experiment objects and correct for GC bias
SE <- SummarizedExperiment(assays = list(counts = counts),
                           rowData = peaks, 
                           colData = DataFrame(names = colnames(counts)))
SE <- addGCBias(SE, genome = BSgenome.Hsapiens.UCSC.hg38)

## Import fine-mapped GWAS posterior probabilities
files <- list.files(path = "CausalDB_Sumstats_Kidney/", full.names = T)
trait.scores <- importBedScore(rowRanges(SE), files, colidx = 5)

## Compute weighted deviation
wDEV_list <- lapply(1:n.runs, function(i) {
  wDEV <- computeWeightedDeviations(SE, trait.scores)
  assays(wDEV)[["z"]]
})

gdev.mat <- apply(abind::abind(wDEV_list, along = 3), c(1,2), function(x) {
  x <- x[!(is.na(x) | is.infinite(x) | is.nan(x))]
  mean(x)
})

## Visualize results (Fig. 6f)
gdev.mat <- t(gdev.mat)
gdev.mat[is.nan(gdev.mat)] <- 0
gdev.mat <- gdev.mat[,apply(gdev.mat, 2, function(x) any(x != 0))]
gdev.mat <- gdev.mat[,c("Gout","Blood_Urea_Nitrogen","Chronic_Kidney_Disease","Chronic_Kidney_Failure","Hypertension","eGFR")]
ggHeat(gdev.mat, heatscale = c(low = 'deepskyblue', mid = 'white', high = '#a50f15'),
       x.lab.size = 10, y.lab.size = 9, dot.highlight.cutoff = qnorm(0.01, lower.tail = F))

##link causal snps to genes 
#Load function
FindLinkedGenes <- function(pheno.cl.use, causal.z, peak.gene.weights) {
  ## Find peaks associated with each trait/cluster pair
  pheno.cluster.peaks <- lapply(pheno.cl.use, function(pheno.cl) {
    p <- strsplit(pheno.cl, split = ": ")[[1]][[1]]
    cl.use <- strsplit(pheno.cl, split = ": ")[[1]][[2]]
    peaks <- unique(as.character(causal.z[[p]][[cl.use]]$peak))
    peaks[peaks %in% rownames(peak.gene.weights)]
  })
  names(pheno.cluster.peaks) <- pheno.cl.use
  
  ## Get coaccessible genes as list
  pheno.cluster.genes <- lapply(pheno.cluster.peaks, function(peaks) {
    peak.gene.weights[peaks,]
  })
  
  
  ## Bind coaccessible gene weights into matrix
  pheno.genes.mat <- do.call(rbind, lapply(pheno.cluster.genes, as.data.frame))
  
  return(pheno.genes.mat)
}

snps <- ImportCausalSNPs(files, colidx = 5)

causal.z <- FindCausalPeaks(SE, trait.scores, min.z = 1, min.PP = 0.025)

pheno.cl.use <- c("Hypertension: M-TAL","Hypertension: C-TAL",
                  "eGFR: M-TAL", "eGFR: C-TAL", 
                  "Chronic_Kidney_Failure: M-TAL", "Chronic_Kidney_Failure: C-TAL")

#Load Cicero connections from "SNARE RNA/AC - Combined Seurat object"
conns <- FilterConns(conns, min.coaccess = 0.1)
conns$Peak1 <- sub("-",":",conns$Peak1)
conns$Peak2 <- sub("-",":",conns$Peak2)

## Use all peaks
unique.peaks <- unique(conns$Peak1)
peaks.gr <- peak2granges(unique.peaks)

peak.gene.weights <- RegionGeneLinks(peaks.gr, conns, link.promoter = T,
                                     promoter.region = c(-3e3, 3e3),
                                     region.name = NULL, weight.col = "coaccess")
peak.gene.weights$region <- sub(":", "-", peak.gene.weights$region)
peak.gene.weights <- as.data.frame(peak.gene.weights %>% group_by(region) %>% top_n(1, weight))
rownames(peak.gene.weights) <- peak.gene.weights$region

#Remove genes not expressed in TAL (snCv3/scCv3)
clusters <- Idents(refquery)
ref <- subset(refquery, cells = rownames(refquery@meta.data[refquery@meta.data$id == "ref",]))
ref.clusters <- Idents(ref)
ref.rna.counts <- GetAssayData(object = ref, slot = "counts")
ref.rna.avg.counts <- GetClusterAvg(ref.rna.counts, ref.clusters, binarize = T)
hist(ref.rna.avg.counts)
quantile(ref.rna.avg.counts)
min.rna.cl.frac <- 0.1
ref.MTAL.expr <- ref.rna.avg.counts[,"M-TAL"]
ref.MTAL.expr <- names(ref.MTAL.expr[ref.MTAL.expr > min.rna.cl.frac])
ref.CTAL.expr <- ref.rna.avg.counts[,"C-TAL"]
ref.CTAL.expr <- names(ref.CTAL.expr[ref.CTAL.expr > min.rna.cl.frac])

query <- subset(refquery, cells = rownames(refquery@meta.data[refquery@meta.data$id == "query",]))
query.clusters <- Idents(query)
query.rna.counts <- GetAssayData(object = query, slot = "counts")
query.rna.avg.counts <- GetClusterAvg(query.rna.counts, query.clusters, binarize = T)
hist(query.rna.avg.counts)
quantile(query.rna.avg.counts)
min.rna.cl.frac <- 0.1
query.MTAL.expr <- query.rna.avg.counts[,"M-TAL"]
query.MTAL.expr <- names(query.MTAL.expr[query.MTAL.expr > min.rna.cl.frac])
query.CTAL.expr <- query.rna.avg.counts[,"C-TAL"]
query.CTAL.expr <- names(query.CTAL.expr[query.CTAL.expr > min.rna.cl.frac])

to.use <- unique(c(ref.MTAL.expr,ref.CTAL.expr,query.MTAL.expr,query.CTAL.expr))
peak.gene.weights <- peak.gene.weights[peak.gene.weights$gene %in% to.use,]
pheno.genes.mat <- FindLinkedGenes(pheno.cl.use, causal.z, peak.gene.weights)

##Find enriched motifs
EM <- FindMotifs(object = KSAC,
                 features = pheno.genes.mat$region)
EM

EM$TF1 <- gsub("\\s*\\([^\\)]+\\)","",EM$motif.name)
EM$TF2 <- gsub(":.*","",EM$TF1)
EM$TF2 <- toupper(EM$TF2)
EM$TF1 <- gsub(".*:","",EM$TF1)
EM$TF1 <- toupper(EM$TF1)

DefaultAssay(KSAC) <- 'ATAC'

#subset to TFs expressed
EM.exp <- EM[EM$TF1 %in% to.use | EM$TF2 %in% to.use,]
MotifPlot(
  object = KSAC,
  motifs = EM.exp$motif[1:5],
  assay = 'ATAC'
)

##Plot TFBS Activities and RNA expression for ESRRB (Fig. 6f)
DefaultAssay(KSAC) <- "chromvar"
KSAC.male <- subset(KSAC, cells = rownames(KSAC@meta.data[KSAC@meta.data$sex == "M",]))
KSAC.female <- subset(KSAC, cells = rownames(KSAC@meta.data[KSAC@meta.data$sex == "F",]))
DotPlot(KSAC.male, features = c("ESRRB"), dot.scale = 6, col.min = 0
) + RotatedAxis() + scale_color_viridis() + DarkTheme()
DotPlot(KSAC.female, features = c("ESRRB"), dot.scale = 6, col.min = 0
) + RotatedAxis() + scale_color_viridis() + DarkTheme()

DefaultAssay(KSAC) <- "RNA"
DefaultAssay(KSAC.female) <- "RNA"
DefaultAssay(KSAC.male) <- "RNA"

DotPlot(KSAC.male, features = c("ESRRB"),
        cols = c("light gray", "dodgerblue4"), dot.scale = 8) + RotatedAxis()
DotPlot(KSAC.female, features = c("ESRRB"),
        cols = c("light gray", "dodgerblue4"), dot.scale = 8) + RotatedAxis()




###TF Network for TAL (ED Fig. 14d)
#Load Cicero connections from "SNARE RNA/AC - Combined Seurat object"
conns.sub <- FilterConns(conns, min.coaccess = 0.1)
conns.sub$Peak1 <- sub("-",":",conns.sub$Peak1)
conns.sub$Peak2 <- sub("-",":",conns.sub$Peak2)

#Identify DEGs from snCv3
KBR.markers <- FindAllMarkers(KBR, only.pos = TRUE, max.cells.per.ident = 500,
                              logfc.threshold = 0.25, min.pct = 0.25)

#Updated function for JASPAR motifs
library(JASPAR2020)
getJasparMotifs <- function(species = "Homo sapiens", 
                            collection = "CORE", ...) {
  opts <- list()
  opts["species"] <- species
  opts["collection"] <- collection
  opts <- c(opts, list(...))
  out <- TFBSTools::getMatrixSet(JASPAR2020::JASPAR2020, opts)
  if (!isTRUE(all.equal(TFBSTools::name(out), names(out)))) 
    names(out) <- paste(names(out), TFBSTools::name(out), sep = "_")
  return(out)
}

#Prepare objects
DefaultAssay(KSAC) <- "ATAC"
Idents(object = KSAC) <- "subclass.l2"
Idents(object = KSAC) <- factor(Idents(object = KSAC), levels = c(
  "POD","dPOD","PEC","PT-S1","PT-S2","PT-S3","aPT","cycPT","dPT","dPT/DTL",
  "DTL2","DTL1","DTL3","dDTL3","ATL","dATL","aTAL1","aTAL2","M-TAL","dM-TAL",
  "C-TAL","dC-TAL","MD","DCT","dDCT","cycDCT","CNT","dCNT","cycCNT","C-PC",
  "M-PC","dM-PC","IMCD","dIMCD","C-IC-A","dC-IC-A","M-IC-A","tPC-IC","IC-B",
  "PapE","EC-GC","EC-AEA","EC-DVR","EC-PTC","dEC-PTC","EC-AVR","dEC","cycEC",
  "EC-LYM","MC","REN","VSMC","VSMC/P","dVSMC","MYOF","cycMYOF","FIB","M-FIB",
  "dM-FIB","aFIB","dFIB","B","PL","T","NKT","MAST","MAC-M2","cycMNP","MDC",
  "cDC","pDC","ncMON","N","SC/NEU"
))

atac.sampled.counts <- GetAssayData(object = KSAC, slot = "counts")
rownames(atac.sampled.counts) <- sub("-",":",rownames(atac.sampled.counts))

min.cells <- 25
atac.sampled.counts <- atac.sampled.counts[rowSums(atac.sampled.counts) > min.cells,]
atac.peaks <- rownames(atac.sampled.counts)
all.peaks <- peak2granges(rownames(atac.sampled.counts))
atac.umap.emb <- Embeddings(object = KSAC, reduction = "umap.atac")
atac.rna.counts <- GetAssayData(object = refquery, slot = "counts")
atac.rna.idents <- Idents(object = refquery)

# Link ATAC peaks to genes (using cicero connections)
peak.gene.df <- RegionGeneLinks(all.peaks, conns.sub, link.promoter = T,
                                promoter.region = c(-3e3, 3e3),
                                region.name = NULL, weight.col = "coaccess")
peak.gene.df <- as.data.frame(peak.gene.df %>% group_by(region) %>% top_n(1, weight))

#Subset to TAL marker genes
TAL.markers <- KBR.markers[KBR.markers$cluster %in% c("C-TAL", "M-TAL"),]$gene
peak.gene.df <- peak.gene.df[peak.gene.df$gene %in% TAL.markers,]
TAL.peaks <- peak2granges(peak.gene.df$region)

# Link ATAC TFs to regions using motifs and correlations
register(MulticoreParam(8))
tf.peak.df <- TFRegionLinks(atac.sampled.counts, atac.umap.emb[colnames(atac.sampled.counts),],
                            TAL.peaks, n.cores = 8)

## Filter TF - region links using a correlation cutoff
hist(tf.peak.df$cor)
quantile(tf.peak.df$cor, na.rm = T)
min.cor <- 0.30
tf.peak.df <- subset(tf.peak.df, cor > min.cor)

## Find TFs expressed in TAL (Expanded to include scCv3)
sn.cells <- rownames(refquery@meta.data[refquery@meta.data$id == "ref",])
sc.cells <- rownames(refquery@meta.data[refquery@meta.data$id == "query",])
sn.rna.avg.counts <- GetClusterAvg(atac.rna.counts[,sn.cells], atac.rna.idents[sn.cells], binarize = T)
sc.rna.avg.counts <- GetClusterAvg(atac.rna.counts[,sc.cells], atac.rna.idents[sc.cells], binarize = T)
hist(sn.rna.avg.counts)
quantile(sn.rna.avg.counts)
min.rna.cl.frac <- 0.2
sn.cTAL.expr <- sn.rna.avg.counts[,"C-TAL"]
sn.mTAL.expr <- sn.rna.avg.counts[,"M-TAL"]
sn.TAL.genes <- c(names(sn.cTAL.expr[sn.cTAL.expr > min.rna.cl.frac]),names(sn.mTAL.expr[sn.mTAL.expr > min.rna.cl.frac]))
sc.cTAL.expr <- sc.rna.avg.counts[,"C-TAL"]
sc.mTAL.expr <- sc.rna.avg.counts[,"M-TAL"]
sc.TAL.genes <- c(names(sc.cTAL.expr[sc.cTAL.expr > min.rna.cl.frac]),names(sc.mTAL.expr[sc.mTAL.expr > min.rna.cl.frac]))
TAL.genes <- unique(c(sn.TAL.genes,sc.TAL.genes))

# Subset to peaks accessible in TAL and TFs expressed in TAL
TAL.networks <- SubsetLinks(tf.peak.df, peak.gene.df, regions = peak.gene.df$region, tfs = TAL.genes)

#Make TF Network
TAL.graph <- graph_from_data_frame(TAL.networks$TF_gene_network)
node.centrality <- eigen_centrality(TAL.graph, directed = F)
node.centrality <- sort(node.centrality$vector, decreasing = T)
ggBarplot(head(node.centrality, n = 20))
tf.gene.df <- TAL.networks$TF_gene_network

ESRRB.neighbors <- subset(tf.gene.df, TF == "ESRRB" | gene == "ESRRB")
ESRRB.network.df <- subset(tf.gene.df, gene %in% ESRRB.neighbors$gene | TF %in% ESRRB.neighbors$TF)
PlotNetwork(ESRRB.network.df, plot.title = "ESRRB Subnetwork", label = T)





###eQTL Enrichment Analyses (ED Fig 14c)
#Load Functions
source("util.func.R") #Trajectory/util.func.R
getEnrichScores <- function(celltype, subclass.l3.markers.sn, subclass.l3.markers.sc, gene.set) {
  subclass.scores <- do.call(rbind, lapply(celltype, function(ct) {
    print(paste("Running for cell type:", ct))
    
    all.de.genes <- unique(c(subclass.l3.markers.sn$gene,subclass.l3.markers.sc$gene))
    sample.de.genes <- unique(c(subclass.l3.markers.sn[subclass.l3.markers.sn$cluster %in% ct,]$gene,
                                subclass.l3.markers.sc[subclass.l3.markers.sc$cluster %in% ct,]$gene))
    subclass.score <- fTest(sample.de.genes, all.de.genes, gene.set)
    subclass.score
  }))
  rownames(subclass.scores) <- celltype
  return(subclass.scores)
}

##Generate conserved marker set for snCv3 and scCv3
#snCv3
subclass.l3.markers.sn <- FindAllMarkers(KBR, features = rownames(KBR),
                                         only.pos = TRUE, max.cells.per.ident = 2000,
                                         logfc.threshold = 0.25, min.pct = 0.1, return.thresh = 0.05)
subclass.l3.markers.sn <- subclass.l3.markers.sn[subclass.l3.markers.sn$p_val_adj < 0.05,]

subclass.l3.markers.sc <- FindAllMarkers(sc.KBR, features = rownames(sc.KBR),
                                         only.pos = TRUE, max.cells.per.ident = 2000,
                                         logfc.threshold = 0.25, min.pct = 0.1, return.thresh = 0.05)
subclass.l3.markers.sc <- subclass.l3.markers.sc[subclass.l3.markers.sc$p_val_adj < 0.05,]


##Calculate enrichment scores for eQTLs
gene.sets.all <- read.delim("eQTL_Gene_Sets_for_Enrichment.txt",sep="\t",header=TRUE)

celltype <- levels(Idents(refquery))
gene.sets <- colnames(gene.sets.all)

gene.set.scores <- do.call(cbind, lapply(gene.sets, function(ct) {
  print(paste("Running for Gene Set:", ct))
  
  gene.set <- gene.sets.all[,ct]
  
  Enrich.scores <- getEnrichScores(celltype, subclass.l3.markers.sn, subclass.l3.markers.sc, gene.set)
  pval <- -log(Enrich.scores$pval,10)
  zval <- scale(pval)
  zval
  
}))

colnames(gene.set.scores) <- gene.sets
rownames(gene.set.scores) <- celltype

ggHeat(gene.set.scores[,1:7], rescaling = "none",clustering = "column",
       heatscale = c(low = '#053061', mid = '#FFFFFF', high = '#67001F'),
       x.lab.size = 10, y.lab.size = 9)
order <- c("Cippa_Acute_Human_Specific","Cippa_Acute_Mouse_Overlap","Cippa_Mid_Acute","Cippa_Late_Human_Specific","Cippa_Late_Mouse_Overlap",
           "Cippa_Late_Fibrosis","Cippa_Late_Recovery")
ggHeat(gene.set.scores[,order], rescaling = "none",clustering = "none",
       heatscale = c(low = '#053061', mid = '#FFFFFF', high = '#67001F'),
       x.lab.size = 10, y.lab.size = 9)



# snCv3 - RNA Velocity (Fig 5a, ED Fig 11b,g-i) ------------------------------------------------------------

###In Python
##Generate loom files from Cellranger outputs. For example:
#velocyto run -m hg38_rmsk.gtf -b KPMP_20191204_10X-R_cellranger_200206/20191204_KPMP_S-1908-000952-R1_SI-GA-H1/outs/filtered_feature_bc_matrix/barcodes.tsv.gz \
#-o Loom_Files/KPMP_20191204A/ \
#KPMP_20191204_10X-R_cellranger_200206/20191204_KPMP_S-1908-000952-R1_SI-GA-H1/outs/possorted_genome_bam.bam \
#refdata-cellranger-GRCh38-3.0.0/genes/genes.gtf


###In R
library(Seurat)
library(velocyto.R)
library(SeuratWrappers)

##Reference sites
#http://htmlpreview.github.io/?https://github.com/satijalab/seurat-wrappers/blob/master/docs/velocity.html
#https://jef.works/blog/2020/08/25/using-scvelo-in-R-using-reticulate/
#http://pklab.med.harvard.edu/velocyto/notebooks/R/SCG71.nb.html

##Prepare Seurat Objects (Correct library IDs for cell barcodes)
#For example:
KB1 <- ReadVelocity(file = "Loom_Files/KPMP_20191204A/possorted_genome_bam_81IDE.loom")
colnames(KB1$spliced) <- gsub("possorted_genome_bam_81IDE:", "", colnames(KB1$spliced))
colnames(KB1$spliced) <- gsub("x", "", colnames(KB1$spliced))
colnames(KB1$unspliced) <- gsub("possorted_genome_bam_81IDE:", "", colnames(KB1$unspliced))
colnames(KB1$unspliced) <- gsub("x", "", colnames(KB1$unspliced))
colnames(KB1$ambiguous) <- gsub("possorted_genome_bam_81IDE:", "", colnames(KB1$ambiguous))
colnames(KB1$ambiguous) <- gsub("x", "", colnames(KB1$ambiguous))
KB1 <- as.Seurat(KB1)


##Merge Seurat objects
KBR.v <- merge(
  x = KB1,
  y = list(KB2,KB3,KB4,KB5,KB6,KB7,KB8,KB9,KB10,KB11,KB12,KB13,KB14,KB15,KB16,
           dKC1,dKC2,dKC3,KC33,KC34,KC35,KC36,KC39,KC40,PA1,KM69,KM71,KM72,KM102,
           KM103,PA3,KC77,KC78,KC79,KBCVD1,KBCVD2,KBCVD3,KBCVD4,KC82,PA4,KC83,KC84,PA5),
  add.cell.ids = c("KB1","KB2","KB3","KB4","KB5","KB6","KB7","KB8","KB9","KB10","KB11","KB12","KB13","KB14","KB15","KB16",
                   "dKC1","dKC2","dKC3","KC33","KC34","KC35","KC36","KC39","KC40","PA1","KM69","KM71","KM72","KM102",
                   "KM103","PA3","KC77","KC78","KC79","KBCVD1","KBCVD2","KBCVD3","KBCVD4","KC82","PA4","KC83","KC84","PA5"),
  merge.data = TRUE
)


##Import the original seurat object
KBR <- LoadH5Seurat("Kidney_Healthy-Injury_Cell_Atlas_snCv3_Seurat_03282022.h5Seurat")
KBR.v <- subset(KBR.v, cells = colnames(KBR))

spliced <- CreateAssayObject(GetAssayData(KBR.v, assay = "spliced"))
unspliced <- CreateAssayObject(GetAssayData(KBR.v, assay = "unspliced"))
ambiguous <- CreateAssayObject(GetAssayData(KBR.v, assay = "ambiguous"))
KBR[["spliced"]] <- spliced
KBR[["unspliced"]] <- unspliced
KBR[["ambiguous"]] <- ambiguous





######Prepare objects for scVelo: aTAL
library(Seurat)
library(velocyto.R)
library(SeuratWrappers)

##Add UMAP and Pseudotime
tal.pseudo <- readRDS("Trajectories/tal.psuedotime.rds")
tal.obj <- readRDS("Trajectories/tal.p2.obj.rds")
length(rownames(tal.pseudo))
cells.use <- rownames(tal.pseudo)[rownames(tal.pseudo) %in% colnames(KBR)]

KBR.TAL <- subset(KBR, cells = cells.use)
umap.coordinates <- tal.obj$embeddings$PCA$umap[cells.use,]
dim(umap.coordinates)

KBR.TAL[["umap.traj"]] <- CreateDimReducObject(embeddings = umap.coordinates, key = "umaptraj_", assay = DefaultAssay(KBR.TAL))

tal.mod <- readRDS("tal.module.assignment.rds")
tal.mod <- tal.mod[rownames(KBR.TAL@meta.data)]
KBR.TAL$modules <- tal.mod

tal.pseudo <- tal.pseudo[rownames(KBR.TAL@meta.data),]
KBR.TAL$pseudotime <- tal.pseudo

DimPlot(KBR.TAL, reduction = "umap.traj", pt.size = 0.5, label = TRUE, group.by = "subclass.l3", repel = TRUE) + NoLegend()
DimPlot(KBR.TAL, reduction = "umap.traj", pt.size = 0.5, label = TRUE, group.by = "clusters", repel = TRUE) + NoLegend()
DimPlot(KBR.TAL, reduction = "umap.traj", pt.size = 0.5, label = TRUE, group.by = "modules", repel = TRUE) + NoLegend()
FeaturePlot(KBR.TAL, features = "pseudotime", reduction = "umap.traj")

##make Anndata object in R
library(reticulate)
use_condaenv("/media/Home_Raid1_Voyager/b1lake/anaconda3/envs/py3", required = TRUE)
scv <- import("scvelo")
scv$logging$print_version()
scv$set_figure_params()

ad <- import("anndata", convert = FALSE)

# pull out spliced and unspliced matrices
expression_matrix <- GetAssayData(KBR.TAL, assay = "RNA", slot = "counts")
expression_matrix <- as.matrix(expression_matrix)
emat <- GetAssayData(KBR.TAL, assay = "spliced", slot = "counts")
emat <- as.matrix(emat)
nmat <- GetAssayData(KBR.TAL, assay = "unspliced", slot = "counts")
nmat <- as.matrix(nmat)
cells <- colnames(emat)
genes <- rownames(expression_matrix)
emat <- emat[genes,]
nmat <- nmat[genes,]

# pull out PCA 
pcs <- Embeddings(KBR.TAL, reduction = "pca")
umap <- Embeddings(KBR.TAL, reduction = "umap.traj")
dfobs <- KBR.TAL@meta.data[,c("clusters","subclass.l3","subclass.l1", "modules", "pseudotime")]
dfvar <- KBR.TAL@assays$RNA@meta.features

adata <- ad$AnnData(
  X=t(expression_matrix),
  obs=dfobs,
  var=dfvar,
  layers=list('spliced'=t(emat), 'unspliced'=t(nmat)),
  obsm=list('X_umap'=umap) 
)

adata$write('scVelo/aTAL_scVelo_annData_object.h5ad', compression='gzip')
#adata = scv$read('scVelo/aTAL_scVelo_annData_object.h5ad')

#See Jupyter Notebook: "snCv3_scVelo_aTAL_Modules.ipynb"



######Prepare objects for scVelo: aPT

##Add UMAP and Pseudotime
#pt.pseudo <- readRDS("Trajectories/pt.psedotime.rds")
pt.obj <- readRDS("Trajectories/pt.p2.obj.rds")
rownames(pt.obj$counts)
#length(rownames(pt.pseudo))
cells.use <- rownames(pt.obj$counts)[rownames(pt.obj$counts) %in% colnames(KBR)]

KBR.pt <- subset(KBR, cells = cells.use)
umap.coordinates <- pt.obj$embeddings$PCA$umap[cells.use,]
dim(umap.coordinates)

KBR.pt[["umap.traj"]] <- CreateDimReducObject(embeddings = umap.coordinates, key = "umaptraj_", assay = DefaultAssay(KBR.pt))

pt.mod <- readRDS("pt.sn.sc.module.assignment.rds")
pt.mod <- pt.mod[rownames(KBR.pt@meta.data)]
KBR.pt$modules <- pt.mod

#pt.pseudo <- pt.pseudo[rownames(KBR.pt@meta.data),]
#KBR.pt$pseudotime <- pt.pseudo

DimPlot(KBR.pt, reduction = "umap.traj", pt.size = 0.5, label = TRUE, group.by = "subclass.l3", repel = TRUE) + NoLegend()
DimPlot(KBR.pt, reduction = "umap.traj", pt.size = 0.5, label = TRUE, group.by = "clusters", repel = TRUE) + NoLegend()
DimPlot(KBR.pt, reduction = "umap.traj", pt.size = 0.5, label = TRUE, group.by = "modules", repel = TRUE) + NoLegend()
#FeaturePlot(KBR.pt, features = "pseudotime", reduction = "umap.traj")

p1 <- DimPlot(KBR.pt, reduction = "umap.traj", pt.size = 0.5, label = TRUE, group.by = "modules",
              repel = TRUE) + NoLegend()
cells.located <- CellSelector(plot = p1)
#"KB1_AGTAACCCAGACTCTA" 


##make Anndata object in R
use_condaenv("/media/Home_Raid1_Voyager/b1lake/anaconda3/envs/py3", required = TRUE)
scv <- import("scvelo")
scv$logging$print_version()
scv$set_figure_params()

ad <- import("anndata", convert = FALSE)

## pull out spliced and unspliced matrices
expression_matrix <- GetAssayData(KBR.pt, assay = "RNA", slot = "counts")
expression_matrix <- as.matrix(expression_matrix)
dim(expression_matrix)
emat <- GetAssayData(KBR.pt, assay = "spliced", slot = "counts")
emat <- as.matrix(emat)
nmat <- GetAssayData(KBR.pt, assay = "unspliced", slot = "counts")
nmat <- as.matrix(nmat)
cells <- colnames(emat)
genes <- rownames(expression_matrix)
emat <- emat[genes,]
dim(emat)
nmat <- nmat[genes,]
dim(nmat)

## pull out PCA 
pcs <- Embeddings(KBR.pt, reduction = "pca")
umap <- Embeddings(KBR.pt, reduction = "umap.traj")
dfobs <- KBR.pt@meta.data[,c("clusters","subclass.l3","subclass.l1", "modules")]
dfvar <- KBR.pt@assays$RNA@meta.features

adata <- ad$AnnData(
  X=t(expression_matrix),
  obs=dfobs,
  var=dfvar,
  layers=list('spliced'=t(emat), 'unspliced'=t(nmat)),
  obsm=list('X_umap'=umap) 
)

adata$write('scVelo/aPT_scVelo_annData_object.h5ad', compression='gzip')
#adata = scv$read('scVelo/aPT_scVelo_annData_object.h5ad')

#See Jupyter Notebook: "snCv3_scVelo_aPT_Modules.ipynb"





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

# snCv3 - Alignment with Mouse AKI Model (ED Fig 7, ED Fig 13f) ----------------------------------------

library(Seurat)
library(SeuratDisk)
library(patchwork)
library(dplyr)
library(Matrix)
library(ggplot2)
library(pagoda2)
library(igraph)

###Prepare Reference Data
KBR #see "snCv3/SNARE2 - Data Integration" 

###Prepare Query Data (Kirita et al., 2020)
hki #see "snCv3 - Marker Genes and Comparison with Reference Data Sets"
hki <- SCTransform(hki, verbose = FALSE)

###Map cell types
#Find anchor genes
anchors <- FindTransferAnchors(
  reference = KBR,
  query = hki,
  normalization.method = "SCT",
  k.filter = NA,
  reference.reduction = "pca",
  reference.neighbors = "pca.annoy.neighbors",
  dims = 1:50
)


hki <- MapQuery(
  anchorset = anchors,
  query = hki,
  reference = KBR,
  refdata = list(
    class = "class",
    subclass.l1 = "subclass.l1",
    subclass.l2 = "subclass.l2",
    subclass.l3 = "subclass.l3"
  ),
  reference.reduction = "pca", 
  reduction.model = "umap"
)


#Subset mouse data to only those with subclass.l1 prediction scores >= 0.8
top.predict <- rownames(hki@meta.data[hki@meta.data$predicted.subclass.l1.score >= 0.8,])
hki <- subset(hki, cells = top.predict)

#Add state annotations
state.map <- KBR@meta.data[!duplicated(KBR@meta.data$subclass.l3),]
rownames(state.map) <- state.map$subclass.l3

state.map <- read.delim("subclass.l3_to_state_mapping.txt")
hki$state.l2 <- hki$predicted.subclass.l3
hki$state.l1 <- hki$predicted.subclass.l3
current.cluster.ids <- state.map$subclass.l3
new.cluster.ids <- state.map$state.l2
hki$state.l2 <- plyr::mapvalues(hki$state.l2, from = current.cluster.ids, to = new.cluster.ids)
new.cluster.ids <- state.map$state.l1
hki$state.l1 <- plyr::mapvalues(hki$state.l1, from = current.cluster.ids, to = new.cluster.ids)




###Plot UMAPs (ED Fig 7a-b)
load("color_factors.robj")

#Subclass level 3
load("color_factors.robj")
Idents(hki) <- "predicted.subclass.l3"
DimPlot(hki, reduction = "ref.umap", label = TRUE, raster=FALSE,
        label.size = 4, repel = TRUE) + ggtitle("Subclass.l3"
        ) + NoLegend() + scale_color_manual(values = alpha(sc.l3.cols[levels(Idents(hki))], 0.4), name = "Subclass"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))

#Original Mouse annotations
load("HKI_color_factors.robj")
Idents(object = hki) <- "celltype"
DimPlot(hki, reduction = "ref.umap", label = TRUE, raster=FALSE,
        label.size = 4, repel = TRUE) + ggtitle("Subclass.l3"
        ) + NoLegend() + scale_color_manual(values = alpha(celltype.cols[levels(Idents(hki))], 0.4), name = "Subclass"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))


#groups
Idents(object = hki) <- "Group"
DimPlot(hki, reduction = "ref.umap", label = FALSE, raster=FALSE,
        label.size = 4, repel = TRUE, split.by = "Group", ncol = 3) + ggtitle("Group"
        ) + scale_color_manual(values = alpha(Group.cols[levels(Idents(object = hki))], 0.4), name = "Group"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))



###State level bar plot (ED Fig 7c)
col.order <- c("Control", "4hours", "12hours", "2days", "14days", "6weeks")
row.order <- c("reference","degenerative","cycling","adaptive - epi","adaptive - str","transitioning")
prop1 <- prop.table(table(hki$state.l2, hki$Group), margin = 2)[row.order,col.order]
barplot(prop1,main = "State Proportions", cex.names = 0.5, 
        names.arg = colnames(prop1),las=2, col = as.character(state.l2.cols[rownames(prop1)]))


###Group level bar plot (ED Fig 7c)
Idents(hki) <- "predicted.subclass.l3"
col.order <- c(
  "POD","dPOD","PEC","PT-S1","PT-S2","PT-S3","aPT","cycPT","dPT","dPT/DTL",
  "DTL2","DTL1","DTL3","dDTL3","ATL","dATL","aTAL1","aTAL2","M-TAL","dM-TAL",
  "C-TAL","dC-TAL","MD","DCT1","DCT2","dDCT","cycDCT","CNT","CNT-PC","dCNT",
  "cycCNT","CCD-PC","OMCD-PC","dOMCD-PC","IMCD","dIMCD","CCD-IC-A","CNT-IC-A",
  "dC-IC-A","OMCD-IC-A","tPC-IC","IC-B","PapE","EC-GC","EC-AEA","EC-DVR","EC-PTC",
  "dEC-PTC","EC-AVR","dEC","cycEC","EC-LYM","MC","REN","VSMC","VSMC/P","dVSMC",
  "MYOF","cycMYOF","FIB","M-FIB","dM-FIB","aFIB","dFIB","B","PL","T","NKT","MAST",
  "MAC-M2","cycMNP","MDC","cDC","pDC","ncMON","N","SC/NEU"
)
col.order <- col.order[col.order %in% levels(Idents(hki))]
Idents(hki) <- factor(Idents(hki), levels = col.order)
col.order <- names(table(Idents(hki))[table(Idents(hki)) > 10])
row.order <- c("Control", "4hours", "12hours", "2days", "14days", "6weeks")
prop1 <- prop.table(table(hki$Group, hki$predicted.subclass.l3), margin = 2)[row.order,col.order]
barplot(prop1,main = "Group Proportions", cex.names = 1, 
        names.arg = colnames(prop1),las=2, col = as.character(Group.cols[rownames(prop1)]))





###Plot marker genes
#TAL genes
modules <- readRDS("mouse.tal.module.assignment.rds")
modules <- modules[rownames(hki@meta.data)]
hki@meta.data$aTAL.modules <- modules

hki$Group.2 <- hki$Group
Idents(hki) <- "Group.2"
hki <- RenameIdents(hki, '12hours' = 'early', '4hours' = 'early', '2days' = 'early',
                    '14days' = 'late', '6weeks' = 'late')
Idents(hki) <- factor(Idents(hki), levels = c("Control","early","late"))
hki$Group.2 <- Idents(hki)

Idents(hki) <- "predicted.subclass.l3"
hki.TAL <- subset(hki, idents = c("C-TAL","aTAL1","aTAL2"))
Idents(hki.TAL) <- "aTAL.modules"
Idents(hki.TAL) <- factor(Idents(hki.TAL), levels = c("blue","yellow","brown","pink"))
hki.TAL <- subset(hki.TAL, idents = c("blue","yellow","brown","pink"))

DotPlot(hki.TAL, features = c("BMP6"), split.by = "Group.2",
        cols = c("#E6ACAC", "#45E6E6", "#804D80"), dot.scale = 8) + RotatedAxis()
DotPlot(hki.TAL, features = c("SLC12A1","EGF","IGF1R"), #scale.max = 30, 
        cols = c("light gray", "dodgerblue4"), dot.scale = 8) + RotatedAxis()
DotPlot(hki.TAL, features = c("SLC12A1","EGF","IGF1R"), split.by = "Group.2",
        cols = c("#E6ACAC", "#45E6E6", "#804D80"), dot.scale = 8) + RotatedAxis()


#STR genes
hki.STR <- subset(hki, idents = c("FIB","aFIB","MYOF"))
Idents(hki.STR) <- factor(Idents(hki.STR), levels = c("FIB","aFIB","MYOF"))
DotPlot(hki.STR, features = c("FBLN5","FLRT2","IGF1","ACVR1","BMPR1A","BMPR1B"), split.by = "Group.2",
        cols = c("#E6ACAC", "#45E6E6", "#804D80"), dot.scale = 8) + RotatedAxis()
DotPlot(hki.STR, features = c("FBLN5","FLRT2","IGF1","ACVR1","BMPR1A","BMPR1B"), 
        cols = c("light gray", "dodgerblue4"), dot.scale = 8) + RotatedAxis()




#Dot plots split by groups.2 excluding control group (ED Fig 13f)
Idents(hki) <- "Group.2"
hki <- subset(hki, idents = c("early","late"))

Idents(hki) <- "predicted.subclass.l3"
hki.TAL <- subset(hki, idents = c("C-TAL","aTAL1","aTAL2"))
Idents(hki.TAL) <- "aTAL.modules"
Idents(hki.TAL) <- factor(Idents(hki.TAL), levels = c("blue","yellow","brown","pink"))
hki.TAL <- subset(hki.TAL, idents = c("blue","yellow","brown","pink"))

DotPlot(hki.TAL, features = c("BMP6"), split.by = "Group.2",
        cols = c("#45E6E6", "#804D80"), dot.scale = 8) + RotatedAxis()
DotPlot(hki.TAL, features = c("SLC12A1","EGF","IGF1R"), split.by = "Group.2",
        cols = c("#45E6E6", "#804D80"), dot.scale = 8) + RotatedAxis()

#STR genes
hki.STR <- subset(hki, idents = c("FIB","aFIB","MYOF"))
Idents(hki.STR) <- factor(Idents(hki.STR), levels = c("FIB","aFIB","MYOF"))
DotPlot(hki.STR, features = c("FBLN5","FLRT2","IGF1","ACVR1","BMPR1A","BMPR1B"), split.by = "Group.2",
        cols = c("#45E6E6", "#804D80"), dot.scale = 8) + RotatedAxis()



# Slide-seq - Figure panels (Fig 2,5, ED Fig 5,13) ------------------------------------------------------
library(ggplot2)
library(Matrix)
library(Giotto)
library(reshape2)
library(patchwork)
library(Seurat)
library(doParallel)
library(foreach)
library(Seurat)
library(dplyr)
library(tidyr)
registerDoParallel(5)

source("GiottoUtils.R")
colorTable <- readRDS("Cluster.Color.Table.RDS")
python_path <- "/anaconda3/envs/giotto_200703/bin/python"

#Giotto Objects 
gioObjFolders <- "/GiottoObjects_norm/"
l1.l2 <- read.table("l1-l2_subclasses.tsv", header = T)

# prepping puck names and addressed
pucks <- foreach(gDir = gioObjFolders) %do% {
  pucks <- list.files(gDir)
  pucks <- pucks[startsWith(pucks, "Puck_") & endsWith(pucks, ".rds")]
  pucks <- substr(pucks, 1, nchar(pucks) - 4)
  pucks
}

# setting up the color table
l1.colors.df <- colorTable[, c('subclass.l1_label', 'subclass.l1_color')]
l1.colors.df <- l1.colors.df[!duplicated(l1.colors.df), ]
l1.colors.df <- l1.colors.df[l1.colors.df$subclass.l1_label %in% l1.l2$name.l1]
l1.colors <- l1.colors.df$subclass.l1_color
names(l1.colors) <- paste0('l1.', l1.colors.df$subclass.l1_label)

l2.colors.df <- colorTable[, c('subclass.l2_label', 'subclass.l2_color')]
l2.colors.df <- l2.colors.df[!duplicated(l2.colors.df), ]
l2.colors.df <- l2.colors.df[l2.colors.df$subclass.l2_label %in% l1.l2$name.l2]
l2.colors <- l2.colors.df$subclass.l2_color
names(l2.colors) <- paste0('l2.', l2.colors.df$subclass.l2_label)

# Giotto's parameters
ginstuc <- createGiottoInstructions(python_path = python_path,
                                    show_plot = F, save_plot = F)

#Load Giotto objects for pucks of interest
pucks.use <- c("Puck_200903_01","Puck_200903_06", "Puck_210119_09", "Puck_210119_18")
ssg.list <- foreach(puck=pucks.use, .final=function(x) setNames(x, pucks.use)) %dopar% {
  i <- if(puck %in% pucks[[1]]) 1 else 2
  ssg <- readRDS(paste0(gioObjFolders[[i]], puck, ".rds"))
  ssg <- subsetGiotto(ssg, cell_ids = ssg@cell_ID[ssg@cell_metadata$has.l1.rctd])
  ssg <- subsetGiotto(ssg, cell_ids = ssg@cell_ID[ssg@cell_metadata$inBoundary])
  ssg 
}







###Plot ED Fig 5c
ssg.puck <- ssg.list[[1]]
ssg.puck <- normalizeGiotto(ssg.puck, scalefactor = 10000, log_norm = T, scale_genes = T, scale_cells = F)
puck = "200903_01"

#Remove outlier beads
spatPlot(ssg.puck)
ssg.puck <- subsetGiottoLocs(ssg.puck, x_min=800, x_max=5700, y_min=300, y_max=5100)

#Plot subclass.l2
sel.cells <- ssg.puck@cell_ID[(ssg.puck@cell_metadata$has.l2.rctd) & 
                                (ssg.puck@cell_metadata$maxWeight.l2 >= 30) & 
                                (substr(ssg.puck@cell_metadata$maxCelltype.l2, 4, 100) %in% l1.l2$name.l2)]
plt3 <- mySpatPlot1(ssg.puck, select_cells = sel.cells, 
                    cell_color = 'maxCelltype.l2', show_legend = T, 
                    cell_color_code = l2.colors, other_cell_color = 'gray95') +
  theme_noAxis(T) 
plt3 <- plt3 + ggtitle(paste0(puck, ' - ', 'maxCelltype.l2')) + 
  theme(plot.title = element_text(face = 'bold', vjust = -3, size = 13))
print(plt3)

#Plot Gloms
to.use <- c("l2.POD","l2.PEC","l2.MC","l2.REN","l2.EC-AEA","l2.EC-GC")
ssg.puck.f <- subsetGiotto(ssg.puck, ssg.puck@cell_ID[ssg.puck@cell_metadata$maxCelltype.l2 %in% to.use])

sel.cells <- ssg.puck.f@cell_ID[(ssg.puck.f@cell_metadata$has.l2.rctd) & 
                                  (ssg.puck.f@cell_metadata$maxWeight.l2 >= 30) & 
                                  (substr(ssg.puck.f@cell_metadata$maxCelltype.l2, 4, 100) %in% l1.l2$name.l2)]
plt5 <- mySpatPlot1(ssg.puck.f, select_cells = sel.cells, 
                    cell_color = 'maxCelltype.l2', show_legend = T, 
                    cell_color_code = l2.colors, other_cell_color = 'gray95') +
  theme_noAxis(T) 
plt5 <- plt5 + ggtitle(paste0(puck, ' - ', 'maxCelltype.l2')) + 
  theme(plot.title = element_text(face = 'bold', vjust = -3, size = 13))
print(plt5)

#Plot tubules
to.use <- c("l2.CNT","l2.M-TAL","l2.M-PC","l2.IMCD","l2.DTL3","l2.DTL1","l2.ATL","l2.DTL2",
            "l2.C-PC","l2.aPT","l2.PT-S3","l2.DCT","l2.aTAL1","l2.C-TAL","l2.MD",
            "l2.PT-S1", "l2.IC-B", "l2.PT-S2","l2.PT-S3","l2.C-IC-A","l2.M-IC-A")
ssg.puck.f <- subsetGiotto(ssg.puck, ssg.puck@cell_ID[ssg.puck@cell_metadata$maxCelltype.l2 %in% to.use])

sel.cells <- ssg.puck.f@cell_ID[(ssg.puck.f@cell_metadata$has.l2.rctd) & 
                                  (ssg.puck.f@cell_metadata$maxWeight.l2 >= 30) & 
                                  (substr(ssg.puck.f@cell_metadata$maxCelltype.l2, 4, 100) %in% l1.l2$name.l2)]
plt5 <- mySpatPlot1(ssg.puck.f, select_cells = sel.cells, 
                    cell_color = 'maxCelltype.l2', show_legend = T, 
                    cell_color_code = l2.colors, other_cell_color = 'gray95') +
  theme_noAxis(T) 
plt5 <- plt5 + ggtitle(paste0(puck, ' - ', 'maxCelltype.l2')) + 
  theme(plot.title = element_text(face = 'bold', vjust = -3, size = 13))
print(plt5)

#Plot Interstitium
to.use <- c("l2.VSMC","l2.VSMC/P","l2.FIB","l2.M-FIB","l2.MYOF","l2.aFIB",
            "l2.EC-PTC","l2.EC-AEA","l2.EC-DVR","l2.EC-AVR","l2.EC-LYM",
            "l2.MAC-M2","l2.T","l2.NKT","l2.N","l2.B","l2.MDC","l2.ncMON","l2.cDC","l2.pDC","l2.MAST","l2.PL")
ssg.puck.f <- subsetGiotto(ssg.puck, ssg.puck@cell_ID[ssg.puck@cell_metadata$maxCelltype.l2 %in% to.use])

sel.cells <- ssg.puck.f@cell_ID[(ssg.puck.f@cell_metadata$has.l2.rctd) & 
                                  (ssg.puck.f@cell_metadata$maxWeight.l2 >= 30) & 
                                  (substr(ssg.puck.f@cell_metadata$maxCelltype.l2, 4, 100) %in% l1.l2$name.l2)]
plt5 <- mySpatPlot1(ssg.puck.f, select_cells = sel.cells, 
                    cell_color = 'maxCelltype.l2', show_legend = T, 
                    cell_color_code = l2.colors, other_cell_color = 'gray95') +
  theme_noAxis(T) 
plt5 <- plt5 + ggtitle(paste0(puck, ' - ', 'maxCelltype.l2')) + 
  theme(plot.title = element_text(face = 'bold', vjust = -3, size = 13))
print(plt5)






###Fibrotic region 1: Plot Fig. 5d and ED Fig. 13a
ssg.puck <- ssg.list[[1]]
ssg.puck <- normalizeGiotto(ssg.puck, scalefactor = 10000, log_norm = T, scale_genes = T, scale_cells = F)
puck = "200903_01"

#Remove outlier beads
spatPlot(ssg.puck)
ssg.puck <- subsetGiottoLocs(ssg.puck, x_min=800, x_max=5700, y_min=300, y_max=5100)

#Plot fibrotic area all subclass.l2 cell types (ED Fig 13a)
spatPlot(ssg.puck)
ssg.puck.f <- subsetGiottoLocs(ssg.puck, x_min=500, x_max=2800, y_min=2300, y_max=4000)

sel.cells <- ssg.puck.f@cell_ID[(ssg.puck.f@cell_metadata$has.l2.rctd) & 
                                  (ssg.puck.f@cell_metadata$maxWeight.l2 >= 30) & 
                                  (substr(ssg.puck.f@cell_metadata$maxCelltype.l2, 4, 100) %in% l1.l2$name.l2)]
plt5 <- mySpatPlot1(ssg.puck.f, select_cells = sel.cells, 
                    cell_color = 'maxCelltype.l2', show_legend = T,
                    point_size = 2,
                    cell_color_code = l2.colors, other_cell_color = 'gray95') +
  theme_noAxis(T) 
plt5 <- plt5 + ggtitle(paste0(puck, ' - ', 'maxCelltype.l2')) + 
  theme(plot.title = element_text(face = 'bold', vjust = -3, size = 13))
print(plt5)

#Plot fibrotic area subclass.l2 main cell types (Fig. 5d)
unique(ssg.puck@cell_metadata$maxCelltype.l2)
to.use <- c("l2.POD","l2.EC-GC","l2.VSMC","l2.FIB","l2.EC-AEA","l2.MYOF","l2.cycMYOF",
            "l2.aFIB", "l2.PL","l2.MAC-M2","l2.MDC","l2.ncMON","l2.cycMNP","l2.NKT","l2.T","l2.cDC",
            "l2.aTAL1","l2.C-TAL")
ssg.puck.f <- subsetGiotto(ssg.puck, ssg.puck@cell_ID[ssg.puck@cell_metadata$maxCelltype.l2 %in% to.use])
ssg.puck.f <- subsetGiottoLocs(ssg.puck.f, x_min=500, x_max=2800, y_min=2300, y_max=4000)

sel.cells <- ssg.puck.f@cell_ID[(ssg.puck.f@cell_metadata$has.l2.rctd) & 
                                  (ssg.puck.f@cell_metadata$maxWeight.l2 >= 30) & 
                                  (substr(ssg.puck.f@cell_metadata$maxCelltype.l2, 4, 100) %in% l1.l2$name.l2)]
plt5 <- mySpatPlot1(ssg.puck.f, select_cells = sel.cells, 
                    cell_color = 'maxCelltype.l2', show_legend = T,
                    point_size = 2,
                    cell_color_code = l2.colors, other_cell_color = 'gray95') +
  theme_noAxis(T) 
plt5 <- plt5 + ggtitle(paste0(puck, ' - ', 'maxCelltype.l2')) + 
  theme(plot.title = element_text(face = 'bold', vjust = -3, size = 13))
print(plt5)

#Plot aTAL genes: ITGB6, EGF, SLC12A1  (Fig. 5d)
plt <- plotGenesRGB_DARK(ssg.puck.f, genes = c('ITGB6', 'EGF', 'SLC12A1'), min.expr = NULL,
                         pt.size = 2,
                         slot='scaled',
                         legend="right")
plt <- plt + ggtitle('ITGB6-EGF-SLC12A1') + theme(plot.title = element_text(face = 'bold', vjust = -3, size = 13))
print(plt)

#Plot MYOF/IMM cell types (Fig. 5d)
plt <- plotMetadata_DARK(ssg.puck.f, metaCols = c("l2.MYOF","l2.aTAL1","l2.T"), 
                         max.expr = NULL, min.expr = c(10,10,10), pt.size = 2,
                         hideAxis = F, legend = "right")
plt <- plt + ggtitle('MYOF-aTAL1-T') + theme(plot.title = element_text(face = 'bold', vjust = -3, size = 13))
print(plt)

#Plot marker genes: CD14, COL1A1, MYH11 (Fig. 5d)
plt <- plotGenesRGB_DARK(ssg.puck.f, genes = c('CD14', 'COL1A1', 'MYH11'), min.expr = NULL,
                         pt.size = 2,
                         slot='scaled',
                         legend="right") 
plt <- plt + ggtitle('CD14-COL1A1-MYH11') + theme(plot.title = element_text(face = 'bold', vjust = -3, size = 13))
print(plt)






###Fibrotic region 2: ED Fig. 13a
ssg.puck <- ssg.list[[3]]
ssg.puck <- normalizeGiotto(ssg.puck, scalefactor = 10000, log_norm = T, scale_genes = T, scale_cells = F)
puck = "210119_09"

#Remove outlier beads
spatPlot(ssg.puck)
ssg.puck <- subsetGiottoLocs(ssg.puck, x_min=500, x_max=5000, y_min=0, y_max=5000)

#Plot fibrotic area for subclass.l2 cell types (ED Fig 13a)
unique(ssg.puck@cell_metadata$maxCelltype.l2)
to.use <- unique(ssg.puck@cell_metadata$maxCelltype.l2)[!unique(ssg.puck@cell_metadata$maxCelltype.l2) %in% "l2.M-TAL"]
ssg.puck.f <- subsetGiotto(ssg.puck, ssg.puck@cell_ID[ssg.puck@cell_metadata$maxCelltype.l2 %in% to.use])
ssg.puck.f <- subsetGiottoLocs(ssg.puck.f, x_min=1500, x_max=3500, y_min=2200, y_max=3500)

plt5 <- mySpatPlot1(ssg.puck.f, select_cells = sel.cells, 
                    cell_color = 'maxCelltype.l2', show_legend = T,
                    point_size = 2,
                    cell_color_code = l2.colors, other_cell_color = 'gray95') +
  theme_noAxis(T) 
plt5 <- plt5 + ggtitle(paste0(puck, ' - ', 'maxCelltype.l2')) + 
  theme(plot.title = element_text(face = 'bold', vjust = -3, size = 13))
print(plt5)






###Fibrotic region 3: ED Fig. 13a
ssg.puck <- ssg.list[[4]]
ssg.puck <- normalizeGiotto(ssg.puck, scalefactor = 10000, log_norm = T, scale_genes = T, scale_cells = F)
puck = "210119_18"

#Remove outlier beads
spatPlot(ssg.puck)
ssg.puck <- subsetGiottoLocs(ssg.puck, x_min=100, x_max=4000, y_min=0, y_max=4100)

#Plot fibrotic area for subclass.l2 cell types (ED Fig 13a)
to.use <- unique(ssg.puck@cell_metadata$maxCelltype.l2)
ssg.puck.f <- subsetGiotto(ssg.puck, ssg.puck@cell_ID[ssg.puck@cell_metadata$maxCelltype.l2 %in% to.use])
ssg.puck.f <- subsetGiottoLocs(ssg.puck.f, x_min=1400, x_max=3500, y_min=3000, y_max=4100)
sel.cells <- ssg.puck.f@cell_ID[(ssg.puck.f@cell_metadata$has.l2.rctd) & 
                                  (ssg.puck.f@cell_metadata$maxWeight.l2 >= 30) & 
                                  (substr(ssg.puck.f@cell_metadata$maxCelltype.l2, 4, 100) %in% l1.l2$name.l2)]
plt5 <- mySpatPlot1(ssg.puck.f, select_cells = sel.cells, 
                    cell_color = 'maxCelltype.l2', show_legend = T,
                    point_size = 2,
                    cell_color_code = l2.colors, other_cell_color = 'gray95') +
  theme_noAxis(T) 
plt5 <- plt5 + ggtitle(paste0(puck, ' - ', 'maxCelltype.l2')) + 
  theme(plot.title = element_text(face = 'bold', vjust = -3, size = 13))
print(plt5)







###Combined dotplot of fibrotic areas (Fig 5d)
#"Puck_200903_01", "Puck_210119_09", "Puck_210119_18"
ssg.puck <- ssg.list[[1]]
ssg.puck.f1 <- subsetGiottoLocs(ssg.puck, x_min=500, x_max=2800, y_min=2300, y_max=4000)

ssg.puck <- ssg.list[[3]]
ssg.puck.f2 <- subsetGiottoLocs(ssg.puck, x_min=1500, x_max=3500, y_min=2200, y_max=3500)

ssg.puck <- ssg.list[[4]]
ssg.puck.f3 <- subsetGiottoLocs(ssg.puck, x_min=1400, x_max=3500, y_min=3000, y_max=4100)

ssg.list.f <- list(ssg.puck.f1,ssg.puck.f2,ssg.puck.f3)


counts.all <- t(merge.sparse(lapply(lapply(ssg.list.f, function(ss) ss@raw_exprs), t)))
meta.all <- do.call(rbind, lapply(ssg.list.f, function(ss) ss@cell_metadata))
ssg.all <- createGiottoObject(raw_exprs = counts.all, cell_metadata = meta.all, instructions = ginstuc)
ssg.all <- subsetGiotto(ssg.all, ssg.all@cell_ID[ssg.all@cell_metadata$inBoundary])

#Plot aTAL markers
genes2plot <- c("ESRRG","EGF","TMEM52B","KCTD1","TFAP2B","ITGB6","BICC1","CST3","TSPAN1","NRP1","TPM1","ERBB2",
                "MET","ARHGAP29","MACC1","PROM1")

rownames(ssg.all@cell_metadata) <- ssg.all@cell_metadata$cell_ID
sss <- CreateSeuratObject(counts = ssg.all@raw_exprs, meta.data = ssg.all@cell_metadata)
sss$maxCelltype.l2 <- substr(sss$maxCelltype.l2, 4, 100)

to.use <- c("aTAL1","C-TAL")
sss <- subset(sss, cells = rownames(sss@meta.data[sss@meta.data$maxCelltype.l2 %in% to.use,]))
sss$maxCelltype.l2 <- factor(sss$maxCelltype.l2, levels = c("aTAL1","C-TAL"))
sss <- SetIdent(sss, value = sss$maxCelltype.l2)

DotPlot(sss, features = (genes2plot), cols = c('gray95', 'orange'), ) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5, size = 9.4))

#Plot aSTR markers
genes2plot <- c("NEGR1","FBLN5","DCN","PDGFRA","LAMA2","FLRT2","LIFR","FGF14","COL16A1","COL14A1",
                "ACTA2","COL1A2","SORBS1",
                "COL1A1")

sss <- CreateSeuratObject(counts = ssg.all@raw_exprs, meta.data = ssg.all@cell_metadata)
sss$maxCelltype.l2 <- substr(sss$maxCelltype.l2, 4, 100)

to.use <- c("FIB","aFIB","MYOF")
sss <- subset(sss, cells = rownames(sss@meta.data[sss@meta.data$maxCelltype.l2 %in% to.use,]))
sss$maxCelltype.l2 <- factor(sss$maxCelltype.l2, levels = c("FIB","aFIB","MYOF"))
sss <- SetIdent(sss, value = sss$maxCelltype.l2)

DotPlot(sss, features = (genes2plot), cols = c('gray95', 'orange'), ) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5, size = 9.4))

#Plot IMM markers
genes2plot <- c("IGKC", "MZB1",
                "IL7R","CD247","CD96", 
                "STAB1", "F13A1", "MRC1", "CD163", 
                "HLA-DRA","CD14","SLC1A3")

sss <- CreateSeuratObject(counts = ssg.all@raw_exprs, meta.data = ssg.all@cell_metadata)
sss$maxCelltype.l2 <- substr(sss$maxCelltype.l2, 4, 100)

to.use <- c("PL","T","MAC-M2","MDC")
sss <- subset(sss, cells = rownames(sss@meta.data[sss@meta.data$maxCelltype.l2 %in% to.use,]))
sss$maxCelltype.l2 <- factor(sss$maxCelltype.l2, levels = c("PL","T","NKT","MAC-M2","MDC","ncMON"))
sss <- SetIdent(sss, value = sss$maxCelltype.l2)

DotPlot(sss, features = (genes2plot), cols = c('gray95', 'orange'), ) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5, size = 9.4))

#Barplot of Imm cell types
order <- c("l2.PL","l2.T","l2.NKT","l2.MAST","l2.MAC-M2","l2.MDC",
           "l2.ncMON","l2.cDC","l2.pDC")
prop1 <- (table(meta.all$maxCelltype.l2))
barplot(prop1[order],main = "Niche Cell Type Counts", cex.names = 0.5, 
        names.arg = order,las=2, col = as.character(l2.colors[order]))




###Full Tissue Puck - 200903_06 (Fig 2c, ED Fig 5c)
ssg.puck <- ssg.list[[2]]
ssg.puck <- normalizeGiotto(ssg.puck, scalefactor = 10000, log_norm = T, scale_genes = T, scale_cells = F)
puck = "200903_06"

#Remove outlier beads
spatPlot(ssg.puck)
ssg.puck <- subsetGiottoLocs(ssg.puck, x_min=800, x_max=6000, y_min=0, y_max=5900)


#Plot subclass.l2 (Fig. 2c, ED Fig 5c)
sel.cells <- ssg.puck@cell_ID[(ssg.puck@cell_metadata$has.l2.rctd) & 
                                (ssg.puck@cell_metadata$maxWeight.l2 >= 30) & 
                                (substr(ssg.puck@cell_metadata$maxCelltype.l2, 4, 100) %in% l1.l2$name.l2)]
plt3 <- mySpatPlot1(ssg.puck, select_cells = sel.cells, 
                    cell_color = 'maxCelltype.l2', show_legend = T, 
                    cell_color_code = l2.colors, other_cell_color = 'gray95') +
  theme_noAxis(T) 
plt3 <- plt3 + ggtitle(paste0(puck, ' - ', 'maxCelltype.l2')) + 
  theme(plot.title = element_text(face = 'bold', vjust = -3, size = 13))
print(plt3)

#Plot Glom (ED Fig 5c)
to.use <- c("l2.POD","l2.PEC","l2.MC","l2.REN","l2.EC-AEA","l2.EC-GC")
ssg.puck.f <- subsetGiotto(ssg.puck, ssg.puck@cell_ID[ssg.puck@cell_metadata$maxCelltype.l2 %in% to.use])

sel.cells <- ssg.puck.f@cell_ID[(ssg.puck.f@cell_metadata$has.l2.rctd) & 
                                  (ssg.puck.f@cell_metadata$maxWeight.l2 >= 30) & 
                                  (substr(ssg.puck.f@cell_metadata$maxCelltype.l2, 4, 100) %in% l1.l2$name.l2)]
plt5 <- mySpatPlot1(ssg.puck.f, select_cells = sel.cells, 
                    cell_color = 'maxCelltype.l2', show_legend = T, 
                    cell_color_code = l2.colors, other_cell_color = 'gray95') +
  theme_noAxis(T) 
plt5 <- plt5 + ggtitle(paste0(puck, ' - ', 'maxCelltype.l2')) + 
  theme(plot.title = element_text(face = 'bold', vjust = -3, size = 13))
print(plt5)

#Plot tubules (ED Fig 5c)
to.use <- c("l2.CNT","l2.M-TAL","l2.M-PC","l2.IMCD","l2.DTL3","l2.DTL1","l2.ATL","l2.DTL2",
            "l2.C-PC","l2.aPT","l2.PT-S3","l2.DCT","l2.aTAL1","l2.C-TAL","l2.MD",
            "l2.PT-S1", "l2.IC-B", "l2.PT-S2","l2.PT-S3","l2.C-IC-A","l2.M-IC-A")
ssg.puck.f <- subsetGiotto(ssg.puck, ssg.puck@cell_ID[ssg.puck@cell_metadata$maxCelltype.l2 %in% to.use])

sel.cells <- ssg.puck.f@cell_ID[(ssg.puck.f@cell_metadata$has.l2.rctd) & 
                                  (ssg.puck.f@cell_metadata$maxWeight.l2 >= 30) & 
                                  (substr(ssg.puck.f@cell_metadata$maxCelltype.l2, 4, 100) %in% l1.l2$name.l2)]
plt5 <- mySpatPlot1(ssg.puck.f, select_cells = sel.cells, 
                    cell_color = 'maxCelltype.l2', show_legend = T, 
                    cell_color_code = l2.colors, other_cell_color = 'gray95') +
  theme_noAxis(T) 
plt5 <- plt5 + ggtitle(paste0(puck, ' - ', 'maxCelltype.l2')) + 
  theme(plot.title = element_text(face = 'bold', vjust = -3, size = 13))
print(plt5)

#Plot Interstitium (ED Fig 5c)
to.use <- c("l2.VSMC","l2.VSMC/P","l2.FIB","l2.M-FIB","l2.MYOF","l2.aFIB",
            "l2.EC-PTC","l2.EC-AEA","l2.EC-DVR","l2.EC-AVR","l2.EC-LYM",
            "l2.MAC-M2","l2.T","l2.NKT","l2.N","l2.B","l2.MDC","l2.ncMON","l2.cDC","l2.pDC","l2.MAST","l2.PL")
ssg.puck.f <- subsetGiotto(ssg.puck, ssg.puck@cell_ID[ssg.puck@cell_metadata$maxCelltype.l2 %in% to.use])

sel.cells <- ssg.puck.f@cell_ID[(ssg.puck.f@cell_metadata$has.l2.rctd) & 
                                  (ssg.puck.f@cell_metadata$maxWeight.l2 >= 30) & 
                                  (substr(ssg.puck.f@cell_metadata$maxCelltype.l2, 4, 100) %in% l1.l2$name.l2)]
plt5 <- mySpatPlot1(ssg.puck.f, select_cells = sel.cells, 
                    cell_color = 'maxCelltype.l2', show_legend = T, 
                    cell_color_code = l2.colors, other_cell_color = 'gray95') +
  theme_noAxis(T) 
plt5 <- plt5 + ggtitle(paste0(puck, ' - ', 'maxCelltype.l2')) + 
  theme(plot.title = element_text(face = 'bold', vjust = -3, size = 13))
print(plt5)



#Plot Renal Corpuscle region (Fig. 2f)
to.use <- c("l2.POD","l2.PEC","l2.MC","l2.REN","l2.MD","l2.EC-AEA","l2.EC-GC")
ssg.puck.f <- subsetGiotto(ssg.puck, ssg.puck@cell_ID[ssg.puck@cell_metadata$maxCelltype.l2 %in% to.use])
ssg.puck.f <- subsetGiottoLocs(ssg.puck.f, x_min=2100, x_max=3200, y_min=3000, y_max=3800)

plt5 <- mySpatPlot1(ssg.puck.f, 
                    cell_color = 'maxCelltype.l2', show_legend = T,
                    point_size = 2,
                    cell_color_code = l2.colors, other_cell_color = 'gray95') +
  theme_noAxis(T) 
plt5 <- plt5 + ggtitle(paste0(puck, ' - ', 'maxCelltype.l2')) + 
  theme(plot.title = element_text(face = 'bold', vjust = -3, size = 13))
print(plt5)

plt <- plotGenesRGB_DARK(ssg.puck.f, genes = c('REN', 'NOS1',"EMCN"), min.expr = NULL,
                         pt.size = 2,
                         slot='scaled',
                         legend="right")
plt <- plt + ggtitle('REN-NOS1-EMCN') + theme(plot.title = element_text(face = 'bold', vjust = -3, size = 13))
print(plt)

#Plot AEA subregion (Fig. 2g)
to.use <- c("l2.POD","l2.EC-GC","l2.VSMC","l2.FIB","l2.EC-AEA","l2.REN",
            "l2.C-TAL","l2.MD","l2.DCT","l2.C-PC",
            "l2.PT-S1","l2.PT-S2","l2.PT-S3")
ssg.puck.f <- subsetGiotto(ssg.puck, ssg.puck@cell_ID[ssg.puck@cell_metadata$maxCelltype.l2 %in% to.use])
ssg.puck.f <- subsetGiottoLocs(ssg.puck.f, x_min=4000, x_max=6000, y_min=2800, y_max=4000)

plt5 <- mySpatPlot1(ssg.puck.f, 
                    cell_color = 'maxCelltype.l2', show_legend = T,
                    point_size = 2,
                    cell_color_code = l2.colors, other_cell_color = 'gray95') +
  theme_noAxis(T) 
plt5 <- plt5 + ggtitle(paste0(puck, ' - ', 'maxCelltype.l2')) + 
  theme(plot.title = element_text(face = 'bold', vjust = -3, size = 13))
print(plt5)

plt <- plotGenesRGB_DARK(ssg.puck.f, genes = c('MYH11', 'SLC5A12',"NPHS2"), min.expr = NULL,
                         pt.size = 2,
                         slot='scaled',
                         legend="right")
plt <- plt + ggtitle('MYH11-SLC5A12-NPHS2') + theme(plot.title = element_text(face = 'bold', vjust = -3, size = 13))
print(plt)






###Plot Medullary pucks (ED Fig 5c)
pucks.use <- c("Puck_210113_37","Puck_210113_28")
ssg.list <- foreach(puck=pucks.use, .final=function(x) setNames(x, pucks.use)) %dopar% {
  i <- if(puck %in% pucks[[1]]) 1 else 2
  ssg <- readRDS(paste0(gioObjFolders[[i]], puck, ".rds"))
  ssg <- subsetGiotto(ssg, cell_ids = ssg@cell_ID[ssg@cell_metadata$has.l1.rctd])
  ssg <- subsetGiotto(ssg, cell_ids = ssg@cell_ID[ssg@cell_metadata$inBoundary])
  ssg 
}

##Outer Medulla - Puck_210113_37
ssg.puck <- ssg.list[[1]]
ssg.puck <- normalizeGiotto(ssg.puck, scalefactor = 10000, log_norm = T, scale_genes = T, scale_cells = F)
puck = "210113_37"

#Remove outlier beads
spatPlot(ssg.puck)
ssg.puck <- subsetGiottoLocs(ssg.puck, x_min=400, x_max=5000, y_min=200, y_max=5000)

#Plot subclass.l2 (ED Fig 5c)
sel.cells <- ssg.puck@cell_ID[(ssg.puck@cell_metadata$has.l2.rctd) & 
                                (ssg.puck@cell_metadata$maxWeight.l2 >= 30) & 
                                (substr(ssg.puck@cell_metadata$maxCelltype.l2, 4, 100) %in% l1.l2$name.l2)]
plt3 <- mySpatPlot1(ssg.puck, select_cells = sel.cells, 
                    cell_color = 'maxCelltype.l2', show_legend = T, 
                    cell_color_code = l2.colors, other_cell_color = 'gray95') +
  theme_noAxis(T) 
plt3 <- plt3 + ggtitle(paste0(puck, ' - ', 'maxCelltype.l2')) + 
  theme(plot.title = element_text(face = 'bold', vjust = -3, size = 13))
print(plt3)

#Plot Gloms (ED Fig 5c)
to.use <- c("l2.POD","l2.PEC","l2.MC","l2.REN","l2.EC-AEA","l2.EC-GC")
ssg.puck.f <- subsetGiotto(ssg.puck, ssg.puck@cell_ID[ssg.puck@cell_metadata$maxCelltype.l2 %in% to.use])

sel.cells <- ssg.puck.f@cell_ID[(ssg.puck.f@cell_metadata$has.l2.rctd) & 
                                  (ssg.puck.f@cell_metadata$maxWeight.l2 >= 30) & 
                                  (substr(ssg.puck.f@cell_metadata$maxCelltype.l2, 4, 100) %in% l1.l2$name.l2)]
plt5 <- mySpatPlot1(ssg.puck.f, select_cells = sel.cells, 
                    cell_color = 'maxCelltype.l2', show_legend = T, 
                    cell_color_code = l2.colors, other_cell_color = 'gray95') +
  theme_noAxis(T) 
plt5 <- plt5 + ggtitle(paste0(puck, ' - ', 'maxCelltype.l2')) + 
  theme(plot.title = element_text(face = 'bold', vjust = -3, size = 13))
print(plt5)

#Plot tubules (ED Fig 5c)
to.use <- c("l2.CNT","l2.M-TAL","l2.M-PC","l2.IMCD","l2.DTL3","l2.DTL1","l2.ATL","l2.DTL2",
            "l2.C-PC","l2.aPT","l2.PT-S3","l2.DCT","l2.aTAL1","l2.C-TAL","l2.MD",
            "l2.PT-S1", "l2.IC-B", "l2.PT-S2","l2.PT-S3","l2.C-IC-A","l2.M-IC-A")
ssg.puck.f <- subsetGiotto(ssg.puck, ssg.puck@cell_ID[ssg.puck@cell_metadata$maxCelltype.l2 %in% to.use])

sel.cells <- ssg.puck.f@cell_ID[(ssg.puck.f@cell_metadata$has.l2.rctd) & 
                                  (ssg.puck.f@cell_metadata$maxWeight.l2 >= 30) & 
                                  (substr(ssg.puck.f@cell_metadata$maxCelltype.l2, 4, 100) %in% l1.l2$name.l2)]
plt5 <- mySpatPlot1(ssg.puck.f, select_cells = sel.cells, 
                    cell_color = 'maxCelltype.l2', show_legend = T, 
                    cell_color_code = l2.colors, other_cell_color = 'gray95') +
  theme_noAxis(T) 
plt5 <- plt5 + ggtitle(paste0(puck, ' - ', 'maxCelltype.l2')) + 
  theme(plot.title = element_text(face = 'bold', vjust = -3, size = 13))
print(plt5)

#Plot Interstitium
to.use <- c("l2.VSMC","l2.VSMC/P","l2.FIB","l2.M-FIB","l2.MYOF","l2.aFIB",
            "l2.EC-PTC","l2.EC-AEA","l2.EC-DVR","l2.EC-AVR","l2.EC-LYM",
            "l2.MAC-M2","l2.T","l2.NKT","l2.N","l2.B","l2.MDC","l2.ncMON","l2.cDC","l2.pDC","l2.MAST","l2.PL")
ssg.puck.f <- subsetGiotto(ssg.puck, ssg.puck@cell_ID[ssg.puck@cell_metadata$maxCelltype.l2 %in% to.use])

sel.cells <- ssg.puck.f@cell_ID[(ssg.puck.f@cell_metadata$has.l2.rctd) & 
                                  (ssg.puck.f@cell_metadata$maxWeight.l2 >= 30) & 
                                  (substr(ssg.puck.f@cell_metadata$maxCelltype.l2, 4, 100) %in% l1.l2$name.l2)]
plt5 <- mySpatPlot1(ssg.puck.f, select_cells = sel.cells, 
                    cell_color = 'maxCelltype.l2', show_legend = T, 
                    cell_color_code = l2.colors, other_cell_color = 'gray95') +
  theme_noAxis(T) 
plt5 <- plt5 + ggtitle(paste0(puck, ' - ', 'maxCelltype.l2')) + 
  theme(plot.title = element_text(face = 'bold', vjust = -3, size = 13))
print(plt5)



##Inner Medulla - Puck_210113_28
ssg.puck <- ssg.list[[2]]
ssg.puck <- normalizeGiotto(ssg.puck, scalefactor = 10000, log_norm = T, scale_genes = T, scale_cells = F)
puck = "210113_28"

#Remove outlier beads
spatPlot(ssg.puck)
ssg.puck <- subsetGiottoLocs(ssg.puck, x_min=300, x_max=5000, y_min=0, y_max=4100)

#Plot subclass.l2 (ED Fig 5c)
sel.cells <- ssg.puck@cell_ID[(ssg.puck@cell_metadata$has.l2.rctd) & 
                                (ssg.puck@cell_metadata$maxWeight.l2 >= 30) & 
                                (substr(ssg.puck@cell_metadata$maxCelltype.l2, 4, 100) %in% l1.l2$name.l2)]
plt3 <- mySpatPlot1(ssg.puck, select_cells = sel.cells, 
                    cell_color = 'maxCelltype.l2', show_legend = T, 
                    cell_color_code = l2.colors, other_cell_color = 'gray95') +
  theme_noAxis(T) 
plt3 <- plt3 + ggtitle(paste0(puck, ' - ', 'maxCelltype.l2')) + 
  theme(plot.title = element_text(face = 'bold', vjust = -3, size = 13))
print(plt3)

#Plot Gloms (ED Fig 5c)
to.use <- c("l2.POD","l2.PEC","l2.MC","l2.REN","l2.EC-AEA","l2.EC-GC")
ssg.puck.f <- subsetGiotto(ssg.puck, ssg.puck@cell_ID[ssg.puck@cell_metadata$maxCelltype.l2 %in% to.use])

sel.cells <- ssg.puck.f@cell_ID[(ssg.puck.f@cell_metadata$has.l2.rctd) & 
                                  (ssg.puck.f@cell_metadata$maxWeight.l2 >= 30) & 
                                  (substr(ssg.puck.f@cell_metadata$maxCelltype.l2, 4, 100) %in% l1.l2$name.l2)]
plt5 <- mySpatPlot1(ssg.puck.f, select_cells = sel.cells, 
                    cell_color = 'maxCelltype.l2', show_legend = T, 
                    cell_color_code = l2.colors, other_cell_color = 'gray95') +
  theme_noAxis(T) 
plt5 <- plt5 + ggtitle(paste0(puck, ' - ', 'maxCelltype.l2')) + 
  theme(plot.title = element_text(face = 'bold', vjust = -3, size = 13))
print(plt5)
 
#Plot tubules (ED Fig 5c)
to.use <- c("l2.CNT","l2.M-TAL","l2.M-PC","l2.IMCD","l2.DTL3","l2.DTL1","l2.ATL","l2.DTL2",
            "l2.C-PC","l2.aPT","l2.PT-S3","l2.DCT","l2.aTAL1","l2.C-TAL","l2.MD",
            "l2.PT-S1", "l2.IC-B", "l2.PT-S2","l2.PT-S3","l2.C-IC-A","l2.M-IC-A")
ssg.puck.f <- subsetGiotto(ssg.puck, ssg.puck@cell_ID[ssg.puck@cell_metadata$maxCelltype.l2 %in% to.use])

sel.cells <- ssg.puck.f@cell_ID[(ssg.puck.f@cell_metadata$has.l2.rctd) & 
                                  (ssg.puck.f@cell_metadata$maxWeight.l2 >= 30) & 
                                  (substr(ssg.puck.f@cell_metadata$maxCelltype.l2, 4, 100) %in% l1.l2$name.l2)]
plt5 <- mySpatPlot1(ssg.puck.f, select_cells = sel.cells, 
                    cell_color = 'maxCelltype.l2', show_legend = T, 
                    cell_color_code = l2.colors, other_cell_color = 'gray95') +
  theme_noAxis(T) 
plt5 <- plt5 + ggtitle(paste0(puck, ' - ', 'maxCelltype.l2')) + 
  theme(plot.title = element_text(face = 'bold', vjust = -3, size = 13))
print(plt5)

#Plot Interstitium (ED Fig 5c)
to.use <- c("l2.VSMC","l2.VSMC/P","l2.FIB","l2.M-FIB","l2.MYOF","l2.aFIB",
            "l2.EC-PTC","l2.EC-AEA","l2.EC-DVR","l2.EC-AVR","l2.EC-LYM",
            "l2.MAC-M2","l2.T","l2.NKT","l2.N","l2.B","l2.MDC","l2.ncMON","l2.cDC","l2.pDC","l2.MAST","l2.PL")
ssg.puck.f <- subsetGiotto(ssg.puck, ssg.puck@cell_ID[ssg.puck@cell_metadata$maxCelltype.l2 %in% to.use])

sel.cells <- ssg.puck.f@cell_ID[(ssg.puck.f@cell_metadata$has.l2.rctd) & 
                                  (ssg.puck.f@cell_metadata$maxWeight.l2 >= 30) & 
                                  (substr(ssg.puck.f@cell_metadata$maxCelltype.l2, 4, 100) %in% l1.l2$name.l2)]
plt5 <- mySpatPlot1(ssg.puck.f, select_cells = sel.cells, 
                    cell_color = 'maxCelltype.l2', show_legend = T, 
                    cell_color_code = l2.colors, other_cell_color = 'gray95') +
  theme_noAxis(T) 
plt5 <- plt5 + ggtitle(paste0(puck, ' - ', 'maxCelltype.l2')) + 
  theme(plot.title = element_text(face = 'bold', vjust = -3, size = 13))
print(plt5)


###




