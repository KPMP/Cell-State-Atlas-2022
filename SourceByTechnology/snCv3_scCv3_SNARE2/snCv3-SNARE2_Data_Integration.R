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

