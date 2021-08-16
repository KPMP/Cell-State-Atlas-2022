## this script is used to generate heatmap and gene modules for figure 5 and 
##  suplemented figure 12

library(Seurat)
library(monocle)
library(VGAM)
library(viridis)
library(dplyr)
library(edgeR)
library(biomaRt)
library(reshape2)
library(slingshot)
library(tradeSeq)
library(nnet)
library(pagoda2)
library(RColorBrewer)
library(WGCNA)
source("./util.func.R")

### module detection
pt.p2 <- readRDS("pt.10X.snare.p2.rds")
sds <- readRDS("pt.sds.3d.rds")
pt.seurat <- get(load("Kidney_10X-RNA_SNARE2-RNA_11302020_Integrated_Seurat_E_PT-only.rda")) 
pt <- as.data.frame(slingPseudotime(sds))

mat <- t(pt.p2$counts)
cellWeights <- slingCurveWeights(sds)
curve1.lambda <- sds@curves$curve1$lambda
curve1.pt <- pt[, 1]
names(curve1.pt) <- rownames(pt)
curve1.pt <- curve1.pt[names(curve1.pt) %in% colnames(mat)]
curve1.pt <- curve1.pt[!is.na(curve1.pt)]

sig.de.genes <- read.table("./obj/pt.sig.de.gene.txt", sep="\t", header=TRUE)
mat <- mat[rownames(mat) %in% sig.de.genes$gene, ]
mat <- mat[, colnames(mat) %in% names(curve1.pt)]

### smoothing
gene_mat <- zoo::rollapply(as.matrix(t(mat)), 30, mean, 
                           align='left',partial=T) %>% apply(2, rescale.and.center, max.quantile=1-1e-3) %>% t
colnames(gene_mat) <- colnames(mat)

### detect modules
powers <- c(c(1:10), seq(from = 12, to=20, by=2));
sft=pickSoftThreshold(t(gene_mat),
                      dataIsExpr = TRUE,powerVector = powers,
                      corFnc = cor, corOptions = list(use = 'p'),networkType = "signed")

sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit, signed R^2",type="n", main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red");

# Red line corresponds to using an R^2 cut-off
abline(h=0.80,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

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

minModuleSize = 20
# Module identification using dynamic tree cut
dynamicMods = cutreeDynamic(dendro = geneTree,
                            method="tree", minClusterSize = minModuleSize)
#dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, method="hybrid", deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = minModuleSize)
#the following command gives the module labels and the size of each module. Lable 0 is reserved for unassigned genes
table(dynamicMods)

dynamicColors = labels2colors(dynamicMods)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", 
                    dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, 
                    guideHang = 0.05, main = "Gene dendrogram and module colors")

merged.moduleColors <- dynamicColors
merged.moduleColors[merged.moduleColors %in% c("red", "pink", "purple")] <- "red"
merged.moduleColors[merged.moduleColors %in% c("blue", "brown", "yellow")] <- "yellow"
 
plotDendroAndColors(geneTree, cbind(dynamicColors, merged.moduleColors),
                    c("original", "merged"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

## extract module gene
module_colors= setdiff(unique(merged.moduleColors), "grey")
all_module_genes <- list()

for(i in 1:length(module_colors)){
  color=module_colors[i]
  module=rownames(gene_mat)[which(merged.moduleColors==color)]
  module.genes <- sig.de.genes[sig.de.genes$gene %in% module, ]
  module.genes$module <- color
  
  all_module_genes[[i]] <- module.genes[order(module.genes$fdr), ]
}
  

## heatmap of top module genes
module.genesets <- lapply(all_module_genes, function(r){
  r
}) %>% dplyr::bind_rows()

### overlay with DEGs
### DEGs based on cluster
pt.p2$getDifferentialGenes(type='PCA', verbose=TRUE, name='cluster', groups=pt.cellannot,
                            upregulated.only=TRUE, append.specificity.metrics=TRUE, append.auc=TRUE)
auc.cutoff <- 0.6
diffgenes <- pt.p2$diffgenes$PCA$cluster
diffgenes <- lapply(1:length(diffgenes), function(r){
  print(r)
  celltype <- names(diffgenes)[r]
  degs <- diffgenes[[r]][diffgenes[[r]]$AUC>auc.cutoff, ]
  if(length(nrow(degs))>0){
    degs$celltype <- celltype
    return(degs)
  }
}) %>% plyr::compact() %>% dplyr::bind_rows()
colnames(diffgenes)[5] <- "gene"

module.genesets$deg <- 0
module.genesets[module.genesets$gene %in% diffgenes$gene, ]$deg <- 1

module.deg <- module.genesets[module.genesets$deg==1, ]
module.deg <- merge(module.deg, diffgenes, by=c("gene"))


gene.list.int <- c("PLSCR1", "DCDC2", "HAVCR1", "AQP1", "VCAM1", "PROM1", "HAVCR1", "LCN2")
#################################################################
# plot example genesets
modules <- unique(module.genesets$module)
ntop <- 10
gene.list <- lapply(1:length(modules), function(r){
  module <- modules[r]
  panel.deg <- module.deg[module.deg$module == module, ]
  panel.deg <- panel.deg[order(panel.deg$AUC, decreasing=T), ]
  panel.deg$gene[1:ntop]
  #module.deg[module.deg$module == module, ]$gene[1:ntop]
}) %>% unlist
gene.list <- gene.list[!is.na(gene.list)]
gene.list <- unique(c(gene.list, "PLSCR1", "DCDC2", "HAVCR1", "AQP1", "VCAM1", "PROM1", "HAVCR1", "LCN2"))


## plotting heatmap
gene_mat <- t(pt.p2$counts)
module.plotsets <- module.genesets[module.genesets$gene %in% gene.list, ]
gene_mat_top <- gene_mat[rownames(gene_mat) %in% module.plotsets$gene, ]

### ordering gene
module.plotsets <- module.plotsets[order(module.plotsets$module), ]

color.order <- c("black", "red", "magenta", "yellow", "turquoise")
module.plotsets$module <- factor(module.plotsets$module, levels=color.order)
module.plotsets <- module.plotsets[order(module.plotsets$module), ]
gene_mat_top <- gene_mat_top[match(module.plotsets$gene, rownames(gene_mat_top)), ]

#gene_mat_top <- gene_mat_top[, order(curve1.pt)]

## ordering cells
curve1.pt <- curve1.pt[order(curve1.pt)]
gene_mat_top <- gene_mat_top[, match(names(curve1.pt), colnames(gene_mat_top))]

## smoothing
gene_mat_top_smooth <- zoo::rollapply(as.matrix(t(gene_mat_top)), 30, mean, align='left',partial=T) %>% apply(2, rescale.and.center, max.quantile=1-1e-3) %>% t
colnames(gene_mat_top_smooth) <- colnames(gene_mat_top)
gene_mat_top <- gene_mat_top_smooth

## heatmap of modules
cellannot <- setNames(snare.tenX.seurat@meta.data$subclass.l3, rownames(snare.tenX.seurat@meta.data))
pt.cellannot <- cellannot[names(cellannot) %in% rownames(snare.tenX.seurat@meta.data)[snare.tenX.seurat@meta.data$subclass.l1 %in% c("PT")]]

color.table <- read.table("./obj/color.table.txt", sep="\t", comment.char="!", header = T)
pt.meta <- snare.tenX.seurat@meta.data
pt.subclass <- data.frame(cells=rownames(pt.meta), subclass=pt.meta$subclass.l2)
subclass.color <- unique(data.frame(subclass=color.table$subclass.l2_label, color=color.table$subclass.l2_color))
pt.subclass.color <- merge(pt.subclass, subclass.color, by=c("subclass"))
colpal <- setNames(pt.subclass.color$color, pt.subclass.color$cells)

cluster.col <- RColorBrewer::brewer.pal(n = 8, name = "Set2")
cluster.col <- data.frame(clusters=seq(1, 8, 1), col=cluster.col)
#pt.cellannot <- factor(setNames(refquery.PT@meta.data$subclass.l3, names(Idents(refquery.PT))), levels=unique(refquery.PT@meta.data$subclass.l3))
pt.cellannot <- pt.cellannot[names(pt.cellannot) %in% colnames(gene_mat_top)]
row.df <- data.frame(clust=as.factor(genesets$module))
col.df <- data.frame(celltype=as.factor(pt.cellannot[colnames(gene_mat_top)]))
rownames(col.df) <- colnames(gene_mat_top)

col.pt <- colpal[names(colpal) %in% names(pt.cellannot)]
col.pt <- col.pt[match(rownames(col.df), names(col.pt))]
names(col.pt) <- col.df[, 1]


col.row <- module.plotsets$module
names(col.row) <- module.plotsets$module


ha <- ComplexHeatmap::HeatmapAnnotation(celltype=col.df[, 1],
                                        which='column', show_annotation_name=FALSE, show_legend=T, border=T, 
                                        col=list(celltype=col.pt))

hl <- ComplexHeatmap::HeatmapAnnotation(cluster=names(col.row),
                                        which='row', show_annotation_name=FALSE, show_legend=T, border=T, 
                                        col=list(cluster=col.row))

hm <- ComplexHeatmap::Heatmap(t(scale(t(gene_mat_top))), cluster_columns = F, cluster_rows = F,show_column_names = F,border=T, 
                              name='expression', show_heatmap_legend = F, show_row_dend = F, 
                              show_row_names=F, top_annotation = ha, left_annotation=hl)

labeled.genes <- rownames(gene_mat_top)[round(seq(1,nrow(gene_mat_top),length.out = 20))]
labeled.genes <- c(labeled.genes, c("PLSCR1","DCDC2","HAVCR1","AQP1","VCAM1","PROM1","HAVCR1","LCN2"))
hm + ComplexHeatmap::rowAnnotation(link = ComplexHeatmap::anno_mark(at = match(labeled.genes, rownames(gene_mat_top)), 
                                                                    labels = labeled.genes, labels_gp = grid::gpar(fontsize = 7)))


## module summary
module.summary <- data.frame(table(module.genesets$module))
names(module.summary) <- c("module", "genes")
module.summary$deg.percent <- 0
for(i in module.summary$module){
  print(i)
  degs <- module.genesets[module.genesets$module==i & module.genesets$deg==1, ]
  module.summary[module.summary$module %in% i, ]$deg.percent <- round(nrow(degs)/nrow(module.genesets[module.genesets$module==i, ]), 2)
}
ggplot(module.summary, aes(x=module, y=genes, fill=module))+
  geom_bar(stat="identity", color="black", alpha=0.7)+
  scale_fill_manual(values=rev(c("yellow", "turquoise", "red", "magenta", "black")))+
  coord_flip() + theme_bw() + theme(legend.position="none") + ylab("# of genes in modules") + xlab("")

### plot mean gene expression of different modules
plists <- list()
color.order <- c("black", "red", "magenta", "yellow", "turquoise")

#par(mfrow=c(length(color.order), 1))
for(i in 1:length(color.order)){
  module <- color.order[i]
  module.gene <- module.plotsets[module.plotsets$module==module, ]$gene
  mean.exp <- data.frame(expression=colMeans(gene_mat_top[rownames(gene_mat_top) %in% module.gene, ]), pseudotime=curve1.pt)
  
  loessMod <- loess(expression ~ pseudotime, data=mean.exp, span=0.50)
  df.curve <- data.frame(pseudotime=mean.exp$pseudotime, fit=loessMod$fitted)
  mean.exp$fit <- df.curve$fit
  
  avg_curve <- list(x=mean.exp$pseudotime, y=mean.exp$fit)
  plot(2,2,xlim=c(0, 12),ylim=c(min(avg_curve$y/max(abs(avg_curve$y))), max(avg_curve$y/max(abs(avg_curve$y)))+0.1),
       type = "n",ylab = " ",xlab=" ",cex.lab=1,yaxt="n",xaxt="n",yaxs="i",xaxs="i")
  
  legend("topright", module, bty="n")
  lines(avg_curve$x, avg_curve$y/max(abs(avg_curve$y)), col=module, lwd=2)
  polygon(c(0,avg_curve$x,14),c(min(avg_curve$y/max(abs(avg_curve$y))),
                                avg_curve$y/max(abs(avg_curve$y)),
                                min(avg_curve$y/max(abs(avg_curve$y)))),
          col=module,border = NA)
}

#### TF motif deviation
snare_atac <- get(load("Kidney_SNARE2_Dual_11302020_Seurat_RNA-AC_Signac-Peaks.rda"))
activity <- snare_atac@assays$Activity@scale.data

## ordering genes
activity.module <- activity[rownames(activity) %in% module.plotsets$gene, ]
module.plotsets.activity <- module.plotsets[module.plotsets$gene %in% rownames(activity.module), ]
activity.module <- activity.module[match(module.plotsets.activity$gene, rownames(activity.module)), ]

## ordering cells
curve1.pt <- curve1.pt[order(curve1.pt)]
curve1.pt.activity <- curve1.pt[names(curve1.pt) %in% colnames(activity.module)]
activity.module.pt <- activity.module[, colnames(activity.module) %in% names(curve1.pt.activity)]
activity.module.pt <- activity.module.pt[, match(names(curve1.pt.activity), colnames(activity.module.pt))]

##plot mean activity
for(i in 1:length(color.order)){
  module <- color.order[i]
  module.gene.activity <- module.plotsets[module.plotsets$module==module, ]$gene
  mean.act <- data.frame(activity=colMeans(activity.module.pt[rownames(activity.module.pt) %in% module.gene.activity, ]), pseudotime=curve1.pt.activity)
  
  loessMod <- loess(activity ~ pseudotime, data=mean.act, span=0.50)
  df.curve <- data.frame(pseudotime=mean.act$pseudotime, fit=loessMod$fitted)
  mean.act$fit <- df.curve$fit
  
  avg_curve <- list(x=mean.act$pseudotime, y=mean.act$fit)
  plot(2,2,xlim=c(0, 12),ylim=c(min(avg_curve$y/max(abs(avg_curve$y))), max(avg_curve$y/max(abs(avg_curve$y)))+0.1),
       type = "n",ylab = " ",xlab=" ",cex.lab=1,yaxt="n",xaxt="n",yaxs="i",xaxs="i")
  
  legend("topright", module, bty="n")
  lines(avg_curve$x, avg_curve$y/max(abs(avg_curve$y)), col=module, lwd=2)
  polygon(c(0,avg_curve$x,14),c(min(avg_curve$y/max(abs(avg_curve$y))),
                                avg_curve$y/max(abs(avg_curve$y)),
                                min(avg_curve$y/max(abs(avg_curve$y)))),
          col=module,border = NA)
}

#### ATAC curves
peaks <- snare_atac@assays$ATAC@counts
pt.peaks <- peaks[, colnames(peaks) %in% colnames(gene_mat_top)]
## order peak matrix
curve1.pt <- curve1.pt[order(curve1.pt)]
curve1.pt.atac <- curve1.pt[names(curve1.pt) %in% colnames(pt.peaks)]
pt.peaks <- pt.peaks[, match(names(curve1.pt.atac), colnames(pt.peaks))]


## assign peak to genes
pt.peaks.annot <- annotPeaks(rownames(pt.peaks))

# plot mean peak signals in each module
for(i in color.order){
  module.gene <- module.genesets[module.genesets$module==i, ]$gene
  module.peaks <- pt.peaks.annot[pt.peaks.annot$SYMBOL %in% module.gene, ]
  module.peaks <- paste(module.peaks$seqnames, module.peaks$start, module.peaks$end, sep="-")
  
  module.peak.mat <- pt.peaks[rownames(pt.peaks) %in% module.peaks, ]
  module.atac.sig.mean <- colMeans(module.peak.mat)
  
  mean.atac <- data.frame(atac.sig=module.atac.sig.mean, pseudotime=curve1.pt.atac)
  
  loessMod <- loess(atac.sig ~ pseudotime, data=mean.atac, span=0.50)
  df.curve <- data.frame(pseudotime=mean.atac$pseudotime, fit=loessMod$fitted)
  mean.atac$fit <- df.curve$fit
  
  avg_curve <- list(x=mean.atac$pseudotime, y=mean.atac$fit)
  plot(2,2,xlim=c(0, 12),ylim=c(min(avg_curve$y/max(abs(avg_curve$y))), max(avg_curve$y/max(abs(avg_curve$y)))+0.1),
       type = "n",ylab = " ",xlab=" ",cex.lab=1,yaxt="n",xaxt="n",yaxs="i",xaxs="i")
  
  legend("topright", i, bty="n")
  lines(avg_curve$x, avg_curve$y/max(abs(avg_curve$y)), col=i, lwd=2)
  polygon(c(0,avg_curve$x,14),c(min(avg_curve$y/max(abs(avg_curve$y))),
                                avg_curve$y/max(abs(avg_curve$y)),
                                min(avg_curve$y/max(abs(avg_curve$y)))),
          col=i,border = NA)
}


####################
### pathways
color.order <- c("black", "red", "magenta", "yellow", "turquoise")

for(i in color.order){
   pathways <- read.csv(paste("./obj/pathways/PT/", i, ".result.curated.csv", sep=""))
  genes <- unique(pathways$Submitted.entities.found)
  
  pathway.unique <- list()
  count <- 1
  for(j in genes){
    pathway <- pathways[pathways$Submitted.entities.found %in% j, ]
    pathway.unique[[count]] <- pathway[1, ]
    #pathway.unique[[count]] <- pathway
    count <- count + 1
  }
  pathway.unique <- dplyr::bind_rows(pathway.unique)
  pathway.unique$logp <- -log(pathway.unique$Entities.pValue)
  pathway.unique <- pathway.unique[order(pathway.unique$logp, decreasing=T), ]
  
  pathway.unique <- pathway.unique[1:nrow(pathway.unique), ]
  print(ggplot(data=pathway.unique, aes(x=reorder(Pathway.name, logp), y=logp)) +
    geom_bar(stat="identity", alpha=0.7, color=i, fill=i) +
    coord_flip() + theme_classic() + theme(legend.position="none") + ylab("- log p-value") + xlab(""))
}

#######################################
## log ratio test for enrichment of fail to recovery genesets
fail.genesets <- read.table("./obj/mouse.fail.to.recovery.sets.txt",
                            sep="\t")

genesets <- c()
for(i in fail.genesets$V1){
  #print(strsplit(i, split=", ")[[1]])
  genesets <- c(genesets, strsplit(i, split=", ")[[1]])
}
genesets <- unique(genesets)

################################
fstat <- lapply(1:length(color.order), function(r){
  module <- color.order[r]
  sample.de.genes <- module.deg[module.deg$module==module, ]$gene
  dev.de.genes <- genesets
  fTest(sample.de.genes, module.genesets$gene, dev.de.genes)
}) %>% dplyr::bind_rows()
fstat$module <- color.order
fstat$module <- factor(color.order, levels=color.order)

ggplot2::ggplot(aes(x=module, y=odds.ratio, fill=module), data=fstat) +
  ggplot2::geom_bar(stat="identity", alpha=0.7) + 
  ggplot2::geom_hline(ggplot2::aes(yintercept=1), color="darkred") +
  ggplot2::scale_y_continuous(limits=c(0, max(fstat$odds.ratio)+0.1), expand=c(0, 0)) +
  ggplot2::labs(x="", y="odd ratio") + theme_bw() + 
  ggplot2::theme(axis.text.x=ggplot2::element_text(angle = 90, hjust = 1, vjust=0.5), legend.position="none",
                 panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank()) +
  scale_fill_manual(values=c("black", "red", "magenta", "yellow", "turquoise"))


##########################################
## extract of module associated cells
gene_mat_top <- readRDS("gene_mat_top.rds")

gene_mat <- pt.p2$counts
module_gene_mat <- gene_mat[, colnames(gene_mat) %in% rownames(gene_mat_top)]
module_gene_mat <- module_gene_mat[rownames(module_gene_mat) %in% colnames(gene_mat_top), ]

reduc <- pcaReduction(t(module_gene_mat), center=TRUE)
reduc.module <- pcaReduction(gene_mat_top, center = TRUE)
knng <- makeKnnGraphReduc(reduc.module, k=30, distance="cosine", center=FALSE, n.cores=10,
                              weight.type='1m')
clusters <- graphCluster(knng)

## remove duplicate if any
reduc <- reduc[-which(rownames(reduc) %in% names(which(duplicated(reduc)))), ]
emb <- getEmb(reduc, method="umap")

# plotting
emb <- getEmb(reduc.module, method="umap")
sccore::embeddingPlot(emb, groups=as.factor(clusters)) + theme_bw()
sccore::embeddingPlot(emb, colors=curve1.pt) + theme_bw()

## plot clusters on original embedding
cellannot <- setNames(snare.tenX.seurat@meta.data$subclass.l3, rownames(snare.tenX.seurat@meta.data))
cellcolors <- setNames(color.table$subclass.l3_color, color.table$subclass.l3_label)
pt.p2$plotEmbedding(type='PCA', embeddingType='umap', groups=cellannot,
                    show.legend=F, mark.clusters=T,shuffle.colors=F,
                    mark.cluster.cex=1, alpha=0.3, group.level.colors =cellcolors,
                    cex=0.2)

pt.p2$plotEmbedding(type='PCA', embeddingType='umap', groups=clusters,
                    show.legend=T, mark.clusters=T,shuffle.colors=F,
                    mark.cluster.cex=1, alpha=0.3,
                    cex=0.2)

sccore::embeddingPlot(pt.p2$embeddings$PCA$umap, colors=curve1.pt, plot.na=F, 
                      size=0.2, alpha=0.3) + theme_bw()

sccore::embeddingPlot(pt.p2$embeddings$PCA$umap, groups=clusters, plot.na=F, 
                      size=0.2, alpha=0.3) + theme_bw()


## check gene expression in each module
gene.colors <- gene_mat_top["PLG", ]
sccore::embeddingPlot(pt.p2$embeddings$PCA$umap, colors=gene.colors, plot.na=F, 
                      size=0.2, alpha=0.3) + theme_bw() + ggtitle("black module gene expression")

gene.colors <- gene_mat_top["MECOM", ]
sccore::embeddingPlot(pt.p2$embeddings$PCA$umap, colors=gene.colors, plot.na=F, 
                      size=0.2, alpha=0.3) + theme_bw() + ggtitle("red module gene expression")
gene.colors <- gene_mat_top["GPX4", ]
sccore::embeddingPlot(pt.p2$embeddings$PCA$umap, colors=gene.colors, plot.na=F, 
                      size=0.2, alpha=0.3) + theme_bw() + ggtitle("red module gene expression")

gene.colors <- gene_mat_top["LRRK2", ]
sccore::embeddingPlot(pt.p2$embeddings$PCA$umap, colors=gene.colors, plot.na=F, 
                      size=0.2, alpha=0.3) + theme_bw() + ggtitle("magenta module gene expression")

### correlation with genesets in modules
## calculate average gene expression in each cluster
cluster.mean.exp <- lapply(unique(clusters), function(r){
  cluster.cells <- names(clusters)[clusters==r]
  cluster.mat <- gene_mat_top[, colnames(gene_mat_top) %in% cluster.cells]
  cluster.mean <- data.frame(rowMeans(cluster.mat))
  colnames(cluster.mean) <- r
  cluster.mean
}) %>% dplyr::bind_cols()

### plot heatmap
annot.row <- data.frame(module=module.plotsets$module)
rownames(annot.row) <- rownames(gene_mat_top)
annotcolors <-  factor(c("black", "red", "magenta", "yellow", "turquoise"),
                                       levels=c("black", "red", "magenta", "yellow", "turquoise"))
names(annotcolors) <- unique(module.plotsets$module)

pheatmap::pheatmap(cluster.mean.exp, cluster_rows=F, annotation_row = annot.row, 
                   annotation_colors=list(module=annotcolors), scale="row", 
                   main="Mean gene expression - cell clusters", border_color = NA)

## get module associated cells
module.cell.clusters <- clusters

module.cell.clusters[module.cell.clusters==9] <- "yellow"
module.cell.clusters[module.cell.clusters==10] <- "turquoise"
module.cell.clusters[module.cell.clusters %in% c(1, 2, 3)] <- "black"
module.cell.clusters[module.cell.clusters %in% c(4)] <- "red"
module.cell.clusters[module.cell.clusters %in% c(5, 6, 7, 8)] <- "magenta"

sccore::embeddingPlot(pt.p2$embeddings$PCA$umap, colors=module.cell.clusters, plot.na=F, 
                      size=0.2, alpha=0.1) + theme_bw() + ggtitle("modules")

### align with single-cell
module.cell.clusters <- readRDS("sn.module.cell.clusters.rds")
kid.integrated.seurat <- get(load("UCSD-WashU-Biopsy_092020_Premiere_Biopsy_Integrated_Seurat_04162021.rda"))

pt.cells <- rownames(kid.integrated.seurat@meta.data[kid.integrated.seurat$subclass.l1 %in% c("PT"), ])
pt.integrated.seurat <- kid.integrated.seurat[, colnames(kid.integrated.seurat) %in% pt.cells]

# remove dPT cycPT cells
removed.cells <- rownames(kid.integrated.seurat@meta.data[kid.integrated.seurat$subclass.l3 %in% c("dPT", "cycPT"), ])
pt.integrated.seurat <- pt.integrated.seurat[, -which(colnames(pt.integrated.seurat) %in% removed.cells)]

pt.integrated.embedding <- pt.integrated.seurat@reductions$ref.umap@cell.embeddings
pt.cellannot <- setNames(pt.integrated.seurat@meta.data$subclass.l3, rownames(pt.integrated.seurat@meta.data))

sc.pt.integrated.embedding <- pt.integrated.embedding[rownames(pt.integrated.embedding) %in% sc.cells, ]
sn.pt.integrated.embedding <- pt.integrated.embedding[rownames(pt.integrated.embedding) %in% sn.cells, ]


clusters <- kmeans(pt.integrated.embedding, centers=20)
clusters <- clusters$cluster

## get knns from integrated embedding
k<- 30
ns <- N2R::Knn(as.matrix(pt.integrated.embedding), k, nThreads=10, indexType='angular')
rownames(ns) <- colnames(ns) <- rownames(pt.integrated.embedding)
sc.cells <- rownames(pt.integrated.seurat@meta.data[pt.integrated.seurat@meta.data$assay %in% c("10X scRNA-seq"), ])
sn.cells <- pt.cells[-which(pt.cells %in% sc.cells)]

module.assignment <- lapply(1:ncol(ns), function(r){
  cell <- colnames(ns)[r]
  if(length(which(cell %in% sn.cells)) >0 ){
    if(length(which(names(module.cell.clusters) %in% cell))>0){
      return(module.cell.clusters[names(module.cell.clusters) %in% cell])
    } else{
      return(NA)
    }
  } else{
    # get neighbor cells
    n.cells <-  names(which(ns[, r] == 1))
    if(length(module.cell.clusters[names(module.cell.clusters) %in% n.cells]) == 0){
      return(NA)
    } else{
      # majority vote
      module.sum <- table(module.cell.clusters[names(module.cell.clusters) %in% n.cells])
      cell.module <- names(module.sum)[which(module.sum == max(module.sum))[1]]
      names(cell.module) <- cell
      return(cell.module)
    }
  }
}) %>% unlist()
module.assignment <- module.assignment[!is.na(module.assignment)]

## plotting check
p1 <- sccore::embeddingPlot(pt.integrated.embedding, groups=module.cell.clusters, plot.na=F,
                      alpha=0.1, size=0.2, colors=module.cell.clusters) + theme_bw() + ggtitle("module assignment")
p2 <- sccore::embeddingPlot(pt.integrated.embedding, colors=module.assignment, plot.na=F,
                      alpha=0.1, size=0.2) + theme_bw()
p3 <- sccore::embeddingPlot(pt.integrated.embedding, groups=as.factor(clusters), plot.na=F,
                            alpha=0.1, size=0.2) + theme_bw()
cowplot::plot_grid(plotlist=list(p1, p2, p3), ncol=2)

### assign cells
sc.cell.assignment <- list()
for(i in 1:20){
  print(i)
  cells <- names(clusters[clusters==i])
  module.sum <- table(module.cell.clusters[names(module.cell.clusters) %in% cells])
  if(length(module.sum)==0){
    sc.cells.cluster <- cells[cells %in% sc.cells]
    sc.cell.assignment[[i]] <- setNames(rep(NA, length(sc.cells.cluster)), sc.cells.cluster)
  } else{
    cell.module <- names(module.sum)[which(module.sum == max(module.sum))[1]]
    sc.cells.cluster <- cells[cells %in% sc.cells]
    sc.cell.assignment[[i]] <- setNames(rep(cell.module, length(sc.cells.cluster)), sc.cells.cluster)
  }
}
sc.cell.assignment <- unlist(sc.cell.assignment)
module.assignment <- c(module.cell.clusters, sc.cell.assignment.removed)

### correlation of modules sc, sn
pt.sc.seurat.integrated <- pt.integrated.seurat[, colnames(pt.integrated.seurat) %in% sc.cells]
pt.sn.seurat.integrated <- pt.integrated.seurat[, colnames(pt.integrated.seurat) %in% sn.cells]

sc.mat <- pt.sc.seurat.integrated@assays$RNA@counts
sc.mat <- sc.mat[rownames(sc.mat) %in% rownames(gene_mat_top), ]
sc.gene_mat <- zoo::rollapply(as.matrix(t(sc.mat)), 30, mean, align='left',partial=T) %>% apply(2, rescale.and.center, max.quantile=1-1e-3) %>% t
colnames(sc.gene_mat) <- colnames(sc.mat)

## remove sc cells
sc.removed <- names(clusters[clusters %in% c(15, 6, 2, 3)])
sc.cell.assignment.removed <- sc.cell.assignment[-which(names(sc.cell.assignment) %in% sc.removed)]

sn.module.mean <- lapply(color.order, function(r){
  cells <- names(module.cell.clusters[module.cell.clusters==r])
  t <- data.frame(rowSums(gene_mat_top[, colnames(gene_mat_top) %in% cells]))
  names(t)[1] <- r
  t
}) %>% dplyr::bind_cols()

sc.module.mean <- lapply(color.order, function(r){
  cells <- names(sc.cell.assignment.removed[sc.cell.assignment.removed==r])
  t <- data.frame(rowSums(sc.gene_mat[, colnames(sc.gene_mat) %in% cells]))
  names(t)[1] <- r
  t
}) %>% dplyr::bind_cols()

sc.module.mean <- sc.module.mean[match(rownames(sn.module.mean), rownames(sc.module.mean)), ]

### align with snare-2
snare.seurat <- get(load("Kidney_10X-RNA_SNARE2-RNA_11302020_Integrated_Seurat_G.rda"))
pt.cells <- rownames(snare.seurat@meta.data[snare.seurat@meta.data$subclass.l2 %in% c("PT-S1", "PT-S2", "PT-S3", "aPT"), ])
snare.seurat.pt <- snare.seurat[, colnames(snare.seurat) %in% pt.cells]
pt.embedding <- snare.seurat.pt@reductions$umap@cell.embeddings

snare.cells <- rownames(snare.seurat.pt@meta.data[snare.seurat.pt@meta.data$assay %in% c("SNARE-Seq2"), ])
pt.embedding.snare <- pt.embedding[rownames(pt.embedding) %in% snare.cells, ]
pt.embedding.tenX <- pt.embedding[-which(rownames(pt.embedding) %in% snare.cells), ]

clusters <- kmeans(pt.embedding, centers=20)
clusters <- clusters$cluster

snare.cell.assignment <- list()
for(i in 1:20){
  print(i)
  cells <- names(clusters[clusters==i])
  module.sum <- table(module.cell.clusters[names(module.cell.clusters) %in% cells])
  if(length(module.sum)==0){
    snare.cells.cluster <- cells[cells %in% snare.cells]
    snare.cell.assignment[[i]] <- setNames(rep(NA, length(snare.cells.cluster)), snare.cells.cluster)
  } else{
    cell.module <- names(module.sum)[which(module.sum == max(module.sum))[1]]
    snare.cells.cluster <- cells[cells %in% snare.cells]
    snare.cell.assignment[[i]] <- setNames(rep(cell.module, length(snare.cells.cluster)), snare.cells.cluster)
  }
}
snare.cell.assignment <- unlist(snare.cell.assignment)
snare.module.assignment <- c(module.cell.clusters, snare.cell.assignment)
saveRDS(snare.module.assignment, "/home/qiwenhu/hubmap/kidney/objects/trajetory/modules/PT/pt.snare.module.assignment.rds")

p1 <- sccore::embeddingPlot(pt.embedding, plot.na=F,
                            alpha=0.1, size=0.2, colors=snare.module.assignment) + theme_bw() + ggtitle("module assignment")
p2 <- sccore::embeddingPlot(pt.embedding, groups=pt.cellannot, plot.na=F,
                            alpha=0.1, size=0.2) + theme_bw()
cowplot::plot_grid(plotlist=list(p1, p2), ncol=2)
