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





