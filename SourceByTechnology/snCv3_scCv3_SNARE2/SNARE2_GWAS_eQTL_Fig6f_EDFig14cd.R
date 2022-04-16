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



