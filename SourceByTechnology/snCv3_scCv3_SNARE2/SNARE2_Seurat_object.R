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





