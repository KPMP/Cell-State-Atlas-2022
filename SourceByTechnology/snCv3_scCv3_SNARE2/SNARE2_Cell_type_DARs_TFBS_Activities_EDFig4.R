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




