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

