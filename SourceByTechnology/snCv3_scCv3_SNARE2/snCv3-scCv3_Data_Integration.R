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



