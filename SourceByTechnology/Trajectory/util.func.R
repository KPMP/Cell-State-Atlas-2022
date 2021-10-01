# related functions for fig 5, supplemental fig 12 and supplemental fig 14
# author: Qiwen Hu - 2021

#' calculate rpkm
rpkm <- function(count, gene_length=NULL){
  total_reads <- colSums(count)
  RPKM <- (count / (gene_length/1000)) %*% diag(1 / (total_reads/1000000))
  colnames(RPKM) <- colnames(count)
  return(RPKM)
}

#' test trajetory associated genes
#' @param condition.vec vector contains condition for each cell
#' @param levels vector of condition levels
test.associated.genes <- function(pt, mat, spline.df=3, n.cores=16, condition.vec=NULL,
                                  levels=NULL){
  require(parallel)
  mat <- mat[, colnames(mat) %in% names(pt)]
  df <- data.frame(do.call(rbind, mclapply(setNames(1:nrow(mat), rownames(mat)), function(i){
    exp <- mat[i,]
    sdf <- data.frame(exp=exp,t=pt[colnames(mat)])
    # model
    m <- mgcv::gam(exp~s(t,k=spline.df),data=sdf,familly=gaussian())
    # background
    m0 <- mgcv::gam(exp~1,data=sdf,familly=gaussian())
    
    fstat <- 0; 
    if(m$deviance>0) {
      fstat <- (deviance(m0) - deviance(m))/(df.residual(m0)-df.residual(m))/(deviance(m)/df.residual(m))
    }
    pval <-  pf(fstat,df.residual(m0)-df.residual(m),df.residual(m),lower.tail = FALSE)
    
    # differences based on conditions
    if(!is.null(condition.vec)){
      sdf <- data.frame(exp=exp,t=pt[colnames(mat)], condition=condition.vec[colnames(mat)])
      #sdf$condition <- factor(sdf$condition, levels=unique(sdf$condition))
      # ordered factor based on condition
      if(is.null(levels)){
        sdf <- mutate(sdf, oCondition = ordered(condition, levels=unique(sdf$condition)))
      } else{
        sdf <- mutate(sdf, oCondition = ordered(condition, levels=levels))
      }
      m.cond <- mgcv::gam(exp ~ oCondition + s(t, k=spline.df) + s(t, by = oCondition), data=sdf,
               family=gaussian())
      stat <- summary(m.cond)$s.table
      return(c(pval=pval, A=diff(range(predict(m))), stat[2:nrow(stat), 4]))
    } else{
      return(c(pval=pval,A=diff(range(predict(m)))))
    }
  }, mc.cores=n.cores,mc.preschedule=T)))
  df$gene <- rownames(mat)
  df$fdr <- p.adjust(df$pval);
  df
}

rescale.and.center <- function(x, center=F, max.quantile=0.99) {
  mx <- quantile(abs(x),max.quantile) # absolute maximum
  if(mx==0) mx<-max(abs(x)) # in case the quantile squashes all the signal
  x[x>mx] <- mx; x[x< -1*mx] <- -1*mx; # trim
  if(center) x <- x-mean(x) # center
  x/max(abs(x)); # scale
}

#' plot gene expression along with psedo-time
#' 
plot_gene_psedotime <- function(mat, gene.name=NULL, pseudotime, groups, dot.size=1, alpha=0.05,
                                clpalatte=NULL, difftest.plot=TRUE, spline.df=3){
  if(is.null(gene.name)){
    stop("please provide gene name to plot")
  }
  
  mat.gene <- mat[gene.name, ]
  groups <- data.frame(cells=names(groups), annot=groups)
  pt <- data.frame(cells=names(pseudotime), pseudotime=pseudotime)
  df <- data.frame(cells=names(mat.gene), exp=mat.gene)
  
  df <- dplyr::left_join(df, groups, by=c("cells"))
  df <- dplyr::left_join(df, pt, by=c("cells"))
  
  # remove NA
  if(length(which(is.na(df$pseudotime)))>0){
    df <- df[!is.na(df$pseudotime), ]
  }
  
  if(difftest.plot){
    df <- mutate(df, oannot = ordered(annot, levels=unique(df$annot)))
    m.cond <- mgcv::gam(exp ~ oannot + s(pseudotime, k=spline.df) + s(pseudotime, by = oannot), data=df,
                        family=gaussian())
    p <- plot(m.cond, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
    diff.df <- data.frame(pseudotime=p[[1]]$x, expression=p[[1]]$fit, condition=levels(df$oannot)[1])
    diff.df <- rbind(diff.df, data.frame(pseudotime=p[[2]]$x, expression=p[[2]]$fit, condition=levels(df$oannot)[2]))
    print(ggplot(diff.df, aes(x=pseudotime, y=expression, group=condition)) + geom_line(aes(color=condition)) +
      theme_bw() + ylab("expression - fitted") + ggtitle(gene.name))
  }
  
  
  if(!is.null(clpalatte)){
    q <- ggplot(aes(pseudotime, exp), data = df) + 
      geom_point(aes_string(color = "annot"), size=dot.size, alpha=alpha, data=df) + 
      scale_color_manual(values=clpalatte)
  } else{
    q <- ggplot(aes(pseudotime, exp), data = df) + 
      geom_point(aes_string(color = "annot"), size=dot.size, alpha=alpha, data=df)
  }
  
  loessMod <- loess(exp ~ pseudotime, data=df, span=0.50)
  df.curve <- data.frame(pseudotime=df$pseudotime, fit=loessMod$fitted)
  df$fit <- df.curve$fit
  
  q <- q + geom_line(aes(x=pseudotime, y=fit), data=df) + guides(colour = guide_legend(override.aes = list(alpha = 1)))+
    theme_bw() + ggtitle(gene.name)
  return(q)
}


#' @param mat count matrix sorted by psedotime, rows are gene, colums are cells
plotPsedotimeHeat <- function (mat, cluster_rows = TRUE, hclust_method = "ward.D2", 
                               num_clusters = 6, hmcols = NULL, add_annotation_row = NULL, 
                               add_annotation_col = NULL, show_rownames = FALSE, 
                               scale_max = 3, scale_min = -3){
  require(pheatmap)
  m <- scale(mat, center=TRUE)
  #m = Matrix::t(scale(Matrix::t(m), center = TRUE))
  #m <- m[is.na(row.names(m)) == FALSE, ]
  m[is.nan(m)] = 0
  m[m > scale_max] = scale_max
  m[m < scale_min] = scale_min
  heatmap_matrix <- m
  
  row_dist <- as.dist((1 - cor(Matrix::t(heatmap_matrix)))/2)
  row_dist[is.na(row_dist)] <- 1
  if (is.null(hmcols)) {
    bks <- seq(-3.1, 3.1, by = 0.1)
    hmcols <- colorRamps::blue2green2red(length(bks) - 1)
  }
  else {
    bks <- seq(-3.1, 3.1, length.out = length(hmcols))
  }
  ph <- pheatmap(heatmap_matrix, useRaster = T, cluster_cols = FALSE, 
                 cluster_rows = cluster_rows, show_rownames = F, show_colnames = F, 
                 clustering_distance_rows = row_dist, clustering_method = hclust_method, 
                 cutree_rows = num_clusters, silent = TRUE, filename = NA, 
                 breaks = bks, border_color = NA, color = hmcols)
  if (cluster_rows) {
    annotation_row <- data.frame(Cluster = factor(cutree(ph$tree_row, 
                                                         num_clusters)))
  }
  else {
    annotation_row <- NULL
  }
  if (!is.null(add_annotation_row)) {
    old_colnames_length <- ncol(annotation_row)
    annotation_row <- cbind(annotation_row, add_annotation_row[row.names(annotation_row), 
    ])
    colnames(annotation_row)[(old_colnames_length + 1):ncol(annotation_row)] <- colnames(add_annotation_row)
  }
  if (!is.null(add_annotation_col)) {
    if (nrow(add_annotation_col) != 100) {
      stop("add_annotation_col should have only 100 rows (check genSmoothCurves before you supply the annotation data)!")
    }
    annotation_col <- add_annotation_col
  }
  else {
    annotation_col <- NA
  }
  
  ph_res <- pheatmap(heatmap_matrix[, ], useRaster = T, cluster_cols = FALSE, 
                     cluster_rows = cluster_rows, show_rownames = show_rownames, 
                     show_colnames = F, clustering_distance_rows = row_dist, 
                     clustering_method = hclust_method, cutree_rows = num_clusters, 
                     annotation_row = annotation_row, annotation_col = annotation_col, 
                     treeheight_row = 20, breaks = bks, fontsize = 6, color = hmcols, 
                     border_color = NA, silent = TRUE, filename = NA)
  grid::grid.rect(gp = grid::gpar("fill", col = NA))
  grid::grid.draw(ph_res$gtable)
  return(ph_res)
  
}

rescale.and.center <- function(x, center=F, max.quantile=0.99) {
  mx <- quantile(abs(x),max.quantile) # absolute maximum
  if(mx==0) mx<-max(abs(x)) # in case the quantile squashes all the signal
  x[x>mx] <- mx; x[x< -1*mx] <- -1*mx; # trim
  if(center) x <- x-mean(x) # center
  x/max(abs(x)); # scale
}

#' clustering on trajectory associated genes
#' @param mat smoothed expression matrix
clmodule <- function(gene_mat, cluster=NULL, nclust=6, recluster=TRUE, ntop=6, col=NULL,
                     plot=TRUE){
  gene_mat.hcl <- hclust(dist(gene_mat), method='ward.D2')  
  gene_mat.hcl.clust <- cutree(gene_mat.hcl, nclust)
  
  clusters <- unique(gene_mat.hcl.clust)
  
  res <- lapply(clusters, function(r){
    genesets <- names(gene_mat.hcl.clust[gene_mat.hcl.clust==r])
    emat <- gene_mat[rownames(gene_mat) %in% genesets, ]
    if(recluster){
      dist <- 1-cor(t(emat))
      hclust.cor <- hclust(as.dist(dist), method="ward.D")
      genesets <- genesets[hclust.cor$order]
    }
    
    genesets.sel <- names(rowSums(emat)[order(rowSums(emat), decreasing=T)])[1:ntop]
    col.df <- data.frame(clust=as.factor(cluster[colnames(gene_mat)]))
    #col.df.colors <- setNames(col[, 2], col[, 1])
    if(!is.null(col)){
      ha <- ComplexHeatmap::HeatmapAnnotation(celltype=col.df[, 1],
                                              which='column', show_annotation_name=FALSE, show_legend=T, border=T, 
                                              col=list(celltype=col))
    } else{
      ha <- ComplexHeatmap::HeatmapAnnotation(type=col.df,
                                              which='column', show_annotation_name=FALSE, show_legend=T, border=T)
    }
    
    hm <- ComplexHeatmap::Heatmap(emat[genesets.sel, ], cluster_columns = F, cluster_rows = F,show_column_names = F,border=T, 
                                name='expression', show_heatmap_legend = F, show_row_dend = F, 
                                show_row_names=F, top_annotation = ha, column_title=paste("mod", r))
    labeled.genes <- genesets.sel
    p <- hm + ComplexHeatmap::rowAnnotation(link = ComplexHeatmap::anno_mark(at = match(genesets.sel, rownames(emat[genesets.sel, ])), 
                                                                        labels = genesets.sel))
    #cowplot::grid.grabExpr(draw(p))
    list(genesets=genesets, emat=emat, plot=p)
  })
  
  ## plotting
  if(plot){
    p <- res[[1]]$plot
    for(i in 2:length(res)){
      p <- p + res[[i]]$plot
    }
   ComplexHeatmap::draw(p)
  }
  return(res)
}

#' do PCA reduction based on given matrix
#' @param mat matrix: rows are genes, column are cells
pcaReduction <- function(mat, nPcs=20, center=FALSE, maxit=1000){
  mat <- t(mat)
  if(center){
    cm <- Matrix::colMeans(mat)
    pcs <- irlba(mat, nv=nPcs, nu=0, center=cm, right_only=FALSE,maxit=maxit,reorth=TRUE)
  } else{
    pcs <- irlba(mat, nv=nPcs, nu=0, right_only=FALSE,maxit=maxit,reorth=TRUE)
  }
  
  if(center) {
    pcs$center <- cm
    pcas <- as.matrix(t(as(t(mat %*% pcs$v), "dgeMatrix") - t(cm %*% pcs$v)))
  } else {
    pcas <- as.matrix(mat %*% pcs$v)
  }
  
  rownames(pcas) <- rownames(mat)
  colnames(pcas) <- paste("PC", seq(1,ncol(pcas)), sep="_")
  return(pcas)
}

#' make knn graph based on reduced spaces
#' @param reduced spaces
makeKnnGraphReduc <- function(reduc, k=30, distance="cosine", center=FALSE, n.cores=10,
                              weight.type='1m'){
  
  if (distance %in% c('cosine','angular')) {
    if (center) {
      reduc <- Matrix::rowMeans(reduc) 
    }
    ns <- N2R::Knn(as.matrix(reduc), k, nThreads=n.cores, indexType='angular')
  } else if (distance == "L2") {
    ns <- N2R::Knn(as.matrix(reduc), k, nThreads=n.cores, indexType='L2')
  } 
  colnames(ns) <- rownames(ns) <- rownames(reduc)
  
  diag(ns) <- 0
  ns <- Matrix::drop0(ns)
  
  # weighting edges - from pagoda2
  x <- reduc
  ns@x <- pmax(0, ns@x)
  if (weight.type=='constant') { 
    ns@x <- 1
  }
  if(weight.type=='1m') { 
    ns@x <- pmax(0, 1-ns@x) 
  }
  
  sns <- (ns+t(ns))/2
  g <- igraph::graph_from_adjacency_matrix(sns, mode='undirected', weighted=TRUE)
  return(g)
}

#' get clusters from graph
graphCluster <- function(g, method="louvain"){
  
  if(method=="louvain"){
    lc <- igraph::cluster_louvain(g)
  }
  
  clusters <- membership(lc)
  comm <- communities(lc)
  
  return(clusters)
}

#' plot cell clusters
#' @param reduc matrix in reduced spaces - rows are cell, column are vector reduced dimentions e.g. PCA
getEmb <- function(reduc, method="umap",
                         dims=2, perplexity=10, n.cores=10){
  if(method == "tSNE"){
    emb <- Rtsne::Rtsne(as.matrix(reduc), dims=dims, is_distance=FALSE, perplexity=perplexity, num_threads=n.cores)$Y
  }
  if(method== "umap"){
    emb <- uwot::umap(as.matrix(reduc), 
                      metric="euclidean", n_threads=n.cores, n_components=dims)
  }
  
  rownames(emb) <- rownames(reduc)
  return(emb)
}


# related functions
randomPatition <- function(cd, k=3){
  # random partition the entire matrix into k part
  cv <- cvTools::cvFolds(nrow(cd), K = k, R = 1)
  cmatrix <- list()
  for(i in 1:k){
    cmatrix[[i]] <- cd[cv$subsets[cv$which == i], ]
  }
  return(cmatrix)
}

# tranform ggdend - modified from cacoa ggend
# dend.data$segments contains col column
ggdendT <- function(dend.data, a = 90) {
  dend.data$segments$colorcode <- paste("colorcode", seq(1, nrow(dend.data$segments)), sep="")
  cols <- setNames(dend.data$segments$col, dend.data$segments$colorcode)
  
  ggplot() +
    geom_segment(data = dend.data$segments, aes(x=x, y=y, xend=xend, yend=yend, col=colorcode), show.legend=FALSE) +
    scale_color_manual(values = cols) + 
    labs(x = "", y = "")  + theme_minimal() +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          panel.grid = element_blank()) +
    geom_text(data = dend.data$labels, aes(x, y, label = label),
              hjust = 1, angle = a, size = 3) + ylim(-0.5, NA)
}

plotPcaSpace <- function(d.counts, d.groups, ref.level, target.level, font.size, palette=NULL){
  bal <- getRndBalances(d.counts)
  pca.res <- prcomp(bal$norm)
  pca.loadings <- bal$psi %*% pca.res$rotation
  
  df.pca <- as.data.frame(pca.res$x)
  
  pc1 <- pca.loadings[,1]
  pc2 <- pca.loadings[,2]
  df.loadings <- as.data.frame(cbind(pc1, pc2) * 10)
  
  
  # ----------- PLOT -----------
  group.names <- c(ref.level, target.level)
  options(repr.plot.width = 15, repr.plot.height = 10)
  rda.plot <- ggplot(df.pca, aes(x=PC1, y=PC2)) +
    #   geom_text(aes(label=rownames(df_pca) %in% samplegroups$trgt),size=4) +
    geom_hline(yintercept=0, linetype="dotted") +
    geom_vline(xintercept=0, linetype="dotted") +
    geom_point(aes(colour = factor(group.names[d.groups + 1] ))) +
    labs(colour="Group") +
    coord_fixed()
  
  if(!is.null(palette)) rda.plot <- rda.plot + scale_fill_manual(values=palette)
  
  dx <- max(df.pca[,'PC1']) - min(df.pca[,'PC1'])
  dy <- max(df.pca[,'PC2']) - min(df.pca[,'PC2'])
  
  
  rda.biplot <- rda.plot +
    geom_segment(data=df.loadings, aes(x=0, xend=pc1, y=0, yend=pc2),
                 color="grey", arrow=arrow(length=unit(0.01,"npc")))  +
    coord_flip(clip = "off") +
    geom_text(data=df.loadings,
              aes(x=pc1,y=pc2,label=rownames(df.loadings)),
              color="black", size=3)
  
  
  
  dx <- max(dx, max(df.loadings[,'pc1']) - min(df.loadings[,'pc1']))
  dy <- max(dy, max(df.loadings[,'pc2']) - min(df.loadings[,'pc2']))
  
  rda.biplot <- rda.biplot + coord_fixed(ratio = dx / dy)
  if(!is.null(font.size)) {
    rda.biplot <- rda.biplot + theme(axis.text=element_text(size=font.size), axis.title=element_text(size=font.size))
  }
  return(rda.biplot)
}

getPcaSpace <- function(d.counts, d.groups, ref.level, target.level){
  bal <- getRndBalances(d.counts)
  pca.res <- prcomp(bal$norm)
  pca.loadings <- bal$psi %*% pca.res$rotation
  df.pca <- as.data.frame(pca.res$x)
  
  pc1 <- pca.loadings[,1]
  pc2 <- pca.loadings[,2]
  df.loadings <- as.data.frame(cbind(pc1, pc2) * 10)
  return(list(pcs=df.pca, df.loadings=df.loadings))
}


classifyCelltype <- function(mat, n.cell.counts=1000, n.iter=1000, p.cutoff=0.001, start.from=1){
  t <- resampleContrast(mat$count, mat$params$d.groups, n.cell.counts=1000, n.iter=1000, start.from=start.from)
  pvals <- getCellSignificance(t$balances)
  sigCelltypes <- names(pvals)[pvals<p.cutoff]
  
  ide <- names(pvals) %in% sigCelltypes
  true <- c(TRUE, rep(FALSE, (mat$params$n.celltypes-1)))
  mcc <- mltools::mcc(ide, true)
  return(list(mcc=mcc, pvals=pvals, sigCelltypes=sigCelltypes))
}

#' extract balance matrix from cacoa object and assign batches in kidney dataset
#' @param samples.region.df data frame contains sampleID, specimen_id, percent.cortex and region
get_balance_matrix <- function(cao, samples.region.df){
  d.counts <- data.frame(anno=cao$cell.groups,
                         group=cao$sample.per.cell[match(names(cao$cell.groups), names(cao$sample.per.cell))]) %>%
    table  %>% rbind %>% t
  d.groups <- cao$sample.groups[rownames(d.counts)] == cao$target.level
  
  
  # get balance matrix
  balance <- getRndBalances(d.counts)$norm %>% as.data.frame()
  sample.region.df <- samples.region.df[match(rownames(balance), samples.region.df[, 1]), ]
  sample.region.df$cat <- "High"
  sample.region.df[sample.region.df$percent.cortex<70, ]$cat <- "Low"
  
  balance <- cbind(data.frame(sample=sample.region.df[, 1]), balance)
  balance$percent.cortex <- sample.region.df$percent.cortex
  balance$region <- sample.region.df$region
  balance$condition <- sample.region.df$condition
  return(balance)
}

#' fit linear models, performance measured by leave-one out
#' @param data data for model fitting - the last column is response
fitLM <- function(data, method="lm"){
  lambda <- 10^seq(-3, 3, length = 100)
  
  predicted <- lapply(1:nrow(data), function(r){
    train <- data[-1*r, ]
    test <- data.frame(data[r, ])
    rownames(test) <- rownames(data)[r]
    if(method == "lm"){
      ml <- lm(percent.cortex ~ ., data=train)
    } else if(method == "lasso"){
      ml <- caret::train(percent.cortex ~., data = train, method = "glmnet",
                         trControl = caret::trainControl("cv", number = 3),
                         tuneGrid=expand.grid(alpha = 1, lambda = lambda)
      )
      
    } else if(method == "ridge"){
      ml <- caret::train(percent.cortex ~., data = train, method = "glmnet",
                         trControl=caret::trainControl("cv", number=3),
                         tuneGrid=expand.grid(alpha=0, lambda=lambda))
      
    } else if(method == "elastic_net"){
      ml <- caret::train(
        percent.cortex ~., data = train, method = "glmnet",
        trControl = caret::trainControl("cv", number = 3),
        tuneLength = 10
      )
    }
    
    if(method == "lm"){
      predict(ml, test[, -ncol(test)])
    } else{
      predict(ml, as.matrix(test[, -ncol(test)]))
    }
  }) %>% unlist()
  
  
  # caculate performance
  SSE <- sum((predicted - data[, ncol(data)])^2)
  SST <- sum((data[, ncol(data)] - mean(data[, ncol(data)]))^2)
  
  rmse = round(sqrt(SSE/nrow(data)), 2)
  corelation <- round(cor(data[, ncol(data)], predicted), 2)
  rsqure <- round(corelation^2, 2)
  
  return(list(cor=corelation, rmse=rmse, rsqure=rsqure, predicted=predicted))
}

#' weighted rank differences
#' rank differenced based on normalized -log(p-values)
#' 
weightedDiff <- function(df1, df2){
  df1$norm <- -log(df1$pvals) / sum(-log(df1$pvals))
  df2$norm <- -log(df2$pvals) / sum(-log(df2$pvals))
  df <- merge(df1, df2, by=c(1))
  df$rankdiff <- df$norm.x - df$norm.y
  
  p <- ggplot(data=df, aes(reorder(celltype, rankdiff), rankdiff)) + geom_bar(stat="identity") + coord_flip() + 
    theme_bw() + ylab("Diff") + xlab("celltype")
  return(p)
}

#' return cda p-values for simulated datasets
cdaPvalues <- function(mat, n.cell.counts=1000, n.iter=1000, p.cutoff=0.001,
                       ref.level="normal", target.level="injury"){
  n.rep <- length(mat$count)
  
  df <- lapply(1:n.rep, function(r){
    t <- resampleContrast(mat$count[[r]], mat$params$d.groups, n.cell.counts=1000, n.iter=1000)
    pvals <- as.data.frame(t(getCellSignificance(t$balances)))
    #sigCelltypes <- names(pvals)[pvals<p.cutoff]
    pvals$rep <- r
    #rownames(pvals) <- ""
    pvals
  }) %>% dplyr::bind_rows()
  
  df$ncell <- mat$params$n.cells
  df$ncelltype <- mat$params$n.celltypes
  df$nsamples <- sum(mat$params$n.samples)
  df$mu <- mat$params$base
  df$mu.pr <- mat$params$clevel
  df$lf <- mat$params$lf
  
  return(df)
}

getBalancePlot <- function(mat, palette, font.size=8, ref.celltype=NULL, p.cutoff=0.1, target.level="disease"){
  bal <- getRndBalances(mat$count)
  #bal <- getRndBalances(mat$count[[1]])
  pca.res <- prcomp(bal$norm)
  t <- resampleContrast(mat$count, mat$params$d.groups, n.cell.counts=1000, n.iter=1000,
                        ref.celltype=ref.celltype)
  #t <- resampleContrast(mat$count[[1]], mat$params$d.groups, n.cell.counts=1000, n.iter=1000,
  #                      ref.celltype=ref.celltype)
  pca.loadings <- bal$psi %*% pca.res$rotation
  
  t$bal <- bal
  t$pca.res <- pca.res
  t$pca.loadings <- pca.loadings
  t$pvals <- getCellSignificance(t$balances)
  
  if(!is.null(ref.celltype)){
    t$pvals[names(t$pvals)==ref.celltype] <- 1
  }
  
  # plot contrast
  p <- plotCellLoadings(t, ordering='by.pvalue', alpha=0.05, signif.threshold=p.cutoff, show.pvals=T,
                        ref.level="Ref", target.level=target.level, font.size=font.size, palette=palette, ref.celltype=ref.celltype)
  return(list(plot=p, cda=t))
}

#' @param count list of count matrices
runConos <- function(count){
  preprocessed <- lapply(count, basicP2proc, n.cores=20, 
                         min.cells.per.gene=0, n.odgenes=2e3, get.largevis=FALSE, make.geneknn=FALSE)
  
  con <- Conos$new(preprocessed, n.cores=20)
  con$buildGraph(k=30, k.self=5, space='PCA', ncomps=30, n.odgenes=2000, 
                 matching.method='mNN', metric='angular', score.component.variance=TRUE, verbose=TRUE)
  con$findCommunities(method=leiden.community, resolution=1)
  #con$embedGraph(method="UMAP", min.dist=0.01, spread=15, n.cores=20, min.prob.lower=1e-3)
  return(con)
}


getBalanceTree <- function(dend, d.counts, d.groups, p.threshold=0.01, ref.level="normal", target.level="injury",
                           colorTree=FALSE){
  log.f <- getLogFreq(d.counts)
  t.tmp <- as.phylo(dend)
  gg.dend <- as.ggdend(dend)
  dend.data <- ggdendro::dendro_data(dend, type = "rectangle")
  dend.data$segments$col <- gg.dend$segments$col
  
  node.pos <- dend.data$segments %$% .[(y == yend) & (yend != 0),]
  node.pos$id <- t.tmp$edge[,1]
  
  innode.pos <- unique(node.pos[,c('x','y','id')])
  
  innode.pos$range <- -1
  for(i in 1:nrow(innode.pos)){
    tmp <- node.pos$xend[node.pos$id == innode.pos$id[i]]
    innode.pos$range[i] <- max(tmp) - min(tmp)
  }
  rownames(innode.pos) <- innode.pos$id
  
  res <- getNodeBalances(t.tmp, log.f)
  balances <- as.data.frame(res$bal)
  colnames(balances) <- rownames(innode.pos)
  
  p.val <- c()
  for(i in 1:ncol(balances)){
    aov.data <- cbind(balances[,i], d.groups)
    colnames(aov.data) <- c('balance', 'group')
    
    res <- aov(balance~group,data=as.data.frame(aov.data))
    p.val <- c(p.val,summary(res)[[1]][1,5])
  }
  # p.val <- p.adjust(p.val, method = 'fdr')
  p.val[is.na(p.val)] <- 1
  
  px <- ggdend(dend.data)
  #px <- as.ggdend(dend)
  
  if(sum(p.val < p.threshold) == 0)
    return(px)
  
  bx <- c()
  by <- c()
  bc <- c()
  
  bx.m <- c()
  by.m <- c()
  bc.m <- c()
  
  pval.x <- c()
  pval.y <- c()
  
  n.range <- c()
  x.range <- c()
  y.range <- c()
  
  idx.target <- d.groups
  for(id.node in 1:ncol(balances)){
    
    if(p.val[id.node] < p.threshold){
      pval.x <- c(pval.x, innode.pos$x[id.node])
      pval.y <- c(pval.y, innode.pos$y[id.node])
    }else
      next
    
    x.tmp <- -balances[,id.node]
    x.tmp <- x.tmp - mean(x.tmp)
    n.range <- c(n.range, max(abs(x.tmp)))  # Ticks
    
    x.tmp <- x.tmp / max(abs(x.tmp)) / 2 * innode.pos$range[id.node] * 0.9
    x.tmp <- x.tmp + innode.pos$x[id.node]
    y.tmp <- idx.target * 0.03 + innode.pos$y[id.node] - 0.05
    
    x.range <- c(x.range, innode.pos$x[id.node] + innode.pos$range[id.node]/2)  # Ticks
    y.range <- c(y.range, innode.pos$y[id.node])  # Ticks
    
    c.tmp <- idx.target
    
    bx <- c(bx, x.tmp)
    by <- c(by, y.tmp)
    bc <- c(bc, c.tmp)
    
    bx.m <- c(bx.m, mean(x.tmp[c.tmp]), mean(x.tmp[!c.tmp]))
    by.m <- c(by.m, mean(y.tmp[c.tmp]), mean(y.tmp[!c.tmp]))
    bc.m <- c(bc.m, TRUE, FALSE)
  }
  
  group.names <- c(ref.level, target.level)
  df.balance.points <- data.frame(x = bx, y = by, group = group.names[bc+1])
  df.balance.mean <- data.frame(x = bx.m, y = by.m, group = group.names[bc.m+1])
  df.balance.range <- data.frame(x = x.range, y = y.range, s = n.range)
  df.pval <- data.frame(x = pval.x, y = pval.y)
  
  
  if(colorTree){
    # add color from dendrogram
    gg.dend <- as.ggdend(dend)
    dend.data$segments$colorcode <- paste("colorcode", seq(1, nrow(dend.data$segments)), sep="")
    cols <- setNames(dend.data$segments$col, dend.data$segments$colorcode)
    cols <- c(cols, "#1F9950", "#CC5252")
    names(cols)[(length(cols)-1):length(cols)] <- c(target.level, ref.level)
    
    px <- px + geom_segment(data = dend.data$segments, aes(x=x, y=y, xend=xend, yend=yend, col=colorcode), show.legend=FALSE) +
      scale_color_manual(values = cols) + geom_point(data = df.balance.points,
                                                     aes(x=x, y=y, col = as.factor(group)), show.legend=FALSE, alpha = 0.5, size = 1) +
      geom_point(data = df.balance.mean,
                 aes(x=x, y=y, col = as.factor(group)),
                 size = 3, shape = 18, stroke = 2, show.legend = FALSE) +
      labs(col=" ") +
      geom_text(data=df.balance.range, mapping=aes(x=x, y=y, label=sprintf('%2.1f',s)), vjust=0, hjust=0)
  } else{
    px <- px + geom_point(data = df.balance.points,
                          aes(x=x, y=y, col = as.factor(group)), alpha = 0.5, size = 1) +
      geom_point(data = df.balance.mean,
                 aes(x=x, y=y, col = as.factor(group)),
                 size = 3, shape = 18, stroke = 2) +
      labs(col=" ") +
      geom_text(data=df.balance.range, mapping=aes(x=x, y=y, label=sprintf('%2.1f',s)), vjust=0, hjust=0)
  }
  
  
  if(length(pval.x) > 0){
    px <- px + geom_point(data = df.pval, aes(x=x, y=y)) +
      # labs(size= sprintf('p-value < %1.2f',p.threshold))
      guides(size = FALSE)
  }
  
  return(px)
}

addInjuryScore <- function(d, meta){
  cbm <- function(d) {
    if(is.leaf(d)){
      nodeCells <- attr(d, "nodesinfo")
      injury.scores <- meta[rownames(meta) %in% nodeCells, ]$injury.score
      
      # caculate mean injury score
      injury.score.mean <- round(mean(injury.scores), 2)
      
      # add attributes in the nodes
      attr(d, "injury.score") <- injury.score.mean
      return(d)
    } else {
      oa <- attributes(d)
      d <- lapply(d,cbm)
      attributes(d) <- oa
      nodeCells <- attr(d, "nodesinfo")
      injury.scores <- meta[rownames(meta) %in% nodeCells, ]$injury.score
      
      injury.score.mean <- round(mean(injury.scores), 2)
      # add attributes in the nodes
      attr(d, "injury.score") <- injury.score.mean
      return(d)
    }
  }
  cbm(d)
}

#' color tree based on the attribute value in the tree
colorDend <- function(dend, attr="clRatio", maxValue=1, 
                      colorpallete=colorRampPalette(c("cyan", "grey", "darkorchid1"))(100)){
  #maxValue <- max(get_nodes_attr(dend, attr))
  #colorpallete=colorRampPalette(c("cyan", "grey", "darkorchid1"))(100)
  
  cc2col <- function(percent, base=0.1){
    cv <- round(percent, 0)
    return(colorpallete[cv])
  }
  
  cbm <- function(dend) {
    if(is.leaf(dend)){
      clRatio <- attr(dend, attr)
      #clPercent <- round(clRatio/maxValue, 2)
      if(clRatio > maxValue){
        clRatio <- 1
      }
      clPercent <- clRatio * 100
      col <- cc2col(clPercent)
      attr(dend, "edgePar") <- c(attr(dend, "edgePar"), list(col=col))
      return(dend)
    } else {
      oa <- attributes(dend)
      dend <- lapply(dend, cbm)
      attributes(dend) <- oa
      clRatio <- attr(dend, attr)
      #clPercent <- round(clRatio/maxValue, 2)
      if(clRatio > maxValue){
        clRatio <- 1
      }
      
      clPercent <- clRatio/maxValue * 100
      col <- cc2col(clPercent)
      attr(dend,"edgePar") <- c(attr(dend,"edgePar"), list(col=col))
      return(dend)
    }
  }
  cbm(dend)
}

#' generate random permuted data from a matrix
permutation <- function(mat, ref=NULL){
  
  #if(is.null(ref)){
  permute.mvalue <- as.vector(mat)
  #} else{
  #  permute.mvalue <- mat[, ref]
  #}
  
  permuted.mat <- lapply(1:nrow(mat), function(r){
    nsample <- ncol(mat)
    sample(permute.mvalue, nsample)
  }) %>% do.call(rbind, .)
  rownames(permuted.mat) <- rownames(mat)
  colnames(permuted.mat) <- colnames(mat)
  
  
  return(permuted.mat)
}



#' @param mat count matrix sorted by psedotime, rows are gene, colums are cells
plotPsedotimeHeat <- function (mat, cluster_rows = TRUE, hclust_method = "ward.D2", 
                               num_clusters = 6, hmcols = NULL, add_annotation_row = NULL, 
                               add_annotation_col = NULL, show_rownames = FALSE, 
                               scale_max = 3, scale_min = -3){
  require(pheatmap)
  m <- scale(mat, center=TRUE)
  m[is.nan(m)] = 0
  m[m > scale_max] = scale_max
  m[m < scale_min] = scale_min
  heatmap_matrix <- m
  
  row_dist <- as.dist((1 - cor(Matrix::t(heatmap_matrix)))/2)
  row_dist[is.na(row_dist)] <- 1
  if (is.null(hmcols)) {
    bks <- seq(-3.1, 3.1, by = 0.1)
    hmcols <- colorRamps::blue2green2red(length(bks) - 1)
  }
  else {
    bks <- seq(-3.1, 3.1, length.out = length(hmcols))
  }
  ph <- pheatmap(heatmap_matrix, useRaster = T, cluster_cols = FALSE, 
                 cluster_rows = cluster_rows, show_rownames = F, show_colnames = F, 
                 clustering_distance_rows = row_dist, clustering_method = hclust_method, 
                 cutree_rows = num_clusters, silent = TRUE, filename = NA, 
                 breaks = bks, border_color = NA, color = hmcols)
  if (cluster_rows) {
    annotation_row <- data.frame(Cluster = factor(cutree(ph$tree_row, 
                                                         num_clusters)))
  }
  else {
    annotation_row <- NULL
  }
  if (!is.null(add_annotation_row)) {
    old_colnames_length <- ncol(annotation_row)
    annotation_row <- cbind(annotation_row, add_annotation_row[row.names(annotation_row), 
    ])
    colnames(annotation_row)[(old_colnames_length + 1):ncol(annotation_row)] <- colnames(add_annotation_row)
  }
  if (!is.null(add_annotation_col)) {
    if (nrow(add_annotation_col) != 100) {
      stop("add_annotation_col should have only 100 rows (check genSmoothCurves before you supply the annotation data)!")
    }
    annotation_col <- add_annotation_col
  }
  else {
    annotation_col <- NA
  }
  
  ph_res <- pheatmap(heatmap_matrix[, ], useRaster = T, cluster_cols = FALSE, 
                     cluster_rows = cluster_rows, show_rownames = show_rownames, 
                     show_colnames = F, clustering_distance_rows = row_dist, 
                     clustering_method = hclust_method, cutree_rows = num_clusters, 
                     annotation_row = annotation_row, annotation_col = annotation_col, 
                     treeheight_row = 20, breaks = bks, fontsize = 6, color = hmcols, 
                     border_color = NA, silent = TRUE, filename = NA)
  grid::grid.rect(gp = grid::gpar("fill", col = NA))
  grid::grid.draw(ph_res$gtable)
  return(ph_res)
  
}

GOanalysis <- function(markers, n, outputFunc="BP"){
  require("org.Hs.eg.db")
  require("Category")
  require("GOstats")
  require("wesanderson")
  
  ENTREZID <- unlist(mget(markers, org.Hs.egSYMBOL2EG, ifnotfound=NA))
  ENTREZID <- ENTREZID[!is.na(ENTREZID)]
  
  
  for(function_type in c("BP", "CC", "MF")){
    
    param <- new("GOHyperGParams", geneIds=ENTREZID,
                 #universe=universe,
                 annotation="org.Hs.eg.db", ontology=function_type,pvalueCutoff=0.1,
                 conditional=FALSE, testDirection="over")
    hyp <- hyperGTest(param)
    sumTable <- summary(hyp)
    
    david=sumTable[1:n,]
    david$p.adj <- p.adjust(david$Pvalue, method="BH")
    david$p.adj=-log(david[,2])
    termNumber=nrow(david)
    david$Term <- factor(david$Term, levels=rev(david$Term))
    
    if(outputFunc==function_type){
      result <- david
    }
    
    p1 <- ggplot(data=david, aes(x=p.adj, y=Term, size=Count, color=p.adj))
    p1 <- p1 + geom_point()
    p1 <- p1 + guides(color = FALSE)
    
    p1 <- p1 + xlab(paste("-log10(Pvalue)", sep="")) + ylab("")
    p1 <- p1 + labs(title=paste("DAVID:", function_type, sep=""))
    p1 <- p1 + theme(axis.text.x=element_text(size=10, face="plain", colour ='black'))
    
    
    p1 <- p1 + theme(axis.text.y=element_text(size=10, face="plain", colour ='black')) + 
      scale_color_gradient(low="blue", high="red")
    
    p1 <- p1 + xlim(min(david$p.adj), max(david$p.adj)) + theme_bw()
    #p1=p1+ylim(-10,+15)
    print(p1)
    #ggsave(file=paste(n,'_' ,function_type, ".png", sep=""), scale=0.8, dpi=600, width = 7, height=1+0.25*termNumber)
  }
  
  return(result)
}


#' sample-level clustering based on different genesets
sample_clustering <- function(rpkm, genesets=NULL, rpkm_cutoff=0.1, num_clusters,
                              patients.df, plot=TRUE, plot.dend=TRUE, num.dim=3,
                              title=NULL){
  sample_expressed <- ncol(rpkm)
  rpkm <- rpkm[(rowSums(rpkm) > rpkm_cutoff*ncol(rpkm)), ]
  rpkm <- rpkm[which(apply(rpkm, 1, min) > 0), ]
  
  # push to monocle object
  sample_sheet <- data.frame(samples=colnames(rpkm))
  rownames(sample_sheet) <- sample_sheet$samples
  gene_annotation <- data.frame(id=rownames(rpkm), gene_short_name=rownames(rpkm))
  rownames(gene_annotation) <- gene_annotation$id
  
  pd <- new("AnnotatedDataFrame", data = sample_sheet)
  fd <- new("AnnotatedDataFrame", data = gene_annotation)
  HSMM <- newCellDataSet(as.matrix(rpkm),
                         phenoData = pd,
                         featureData = fd,
                         lowerDetectionLimit = 0.1,
                         expressionFamily = tobit(Lower = 0.1))
  
  HSMM <- estimateSizeFactors(HSMM)
  HSMM <- detectGenes(HSMM, min_expr = 0.1)
  
  if(!is.null(genesets)){
    expressed_genes <- row.names(subset(fData(HSMM),
                                        num_cells_expressed >= 1))
  } else{
    expressed_genes <- row.names(subset(fData(HSMM),
                                        num_cells_expressed == nrow(sample_sheet)))
  }
  
  # reduce dimension and plot cells
  HSMM <- reduceDimension(HSMM[expressed_genes, ], max_components = 2, num_dim = num.dim,
                          reduction_method = 'tSNE', verbose = TRUE, perplexity= 3)
  HSMM <- clusterCells(HSMM, num_clusters=num_clusters)
  
  conditions <- setNames(patients.df$condition.l1, patients.df$patient_id)
  pData(HSMM)$condition.l2 <- patients.df$condition.l2
  pData(HSMM)$tissue_type <- patients.df$tissue_type
  pData(HSMM)$sex <- patients.df$sex
  pData(HSMM)$condition <- conditions
  
  if(plot){
    p1 <- plot_cell_clusters(HSMM)
    p2 <- plot_cell_clusters(HSMM, color_by="condition")
    p3 <- plot_cell_clusters(HSMM, color_by="condition.l2")
    p4 <- plot_cell_clusters(HSMM, color_by="tissue_type", show_cell_names=TRUE, cell_name_size=1)
    p5 <- plot_cell_clusters(HSMM, color_by="sex", show_cell_names=TRUE, cell_name_size=1)
    plotlist <- list(p1, p2, p3, p4, p5)
  } else{
    plotlist=NULL
  }
  
  if(plot.dend){
    sample.cor <- cor(HSMM[expressed_genes, ]@assayData$exprs)
    dend <- hclust(as.dist(sample.cor), method='ward.D2')
    dend <- as.dendrogram(dend)
    
    leaf.labels <- labels(dend)
    leaf.condition <- patients.df[match(leaf.labels, patients.df$specimen), ]$condition.l1
    leaf.condition[leaf.condition == "AKI"] <- "red"
    leaf.condition[leaf.condition == "Ref"] <- "blue"
    leaf.condition[leaf.condition == "CKD"] <- "green"
    
    print(dend %>% set("leaves_pch", 19) %>%  # node point type
            set("leaves_cex", 1) %>%  # node point size
            set("leaves_col", leaf.condition) %>% # node point color
            set("labels_cex", 0.7) %>% hang.dendrogram(hang=0.5) %>%
            plot(main = title))
  }
  
  return(list(monocle=HSMM, plots=plotlist))
}

quickUmap <- function(reduction, space="PCA", dims=NULL, n.neighbors=30L, n.components=2L,
                      metric='cosine', n.epochs=NULL, learning.rate=1.0, min.dist=0.3,
                      spread=1.0, set.op.mix.ratio=1.0, local.connectivity=1L, repulsion.strength=1,
                      negative.sample.rate=5, a=NULL, b=NULL, uwot.sgd=FALSE, seed.use=22,
                      metric.kwds=NULL, angular.rp.forest=FALSE, verbose=TRUE){
  require(uwot)
  if(is.null(reduction)){
    stop("please run dimention reduction first")
  }
  
  if(!is.null(dims)){
    if(max(dims) > ncol(reduction)){
      warning("number of dims is larger than current variables, setting dims equal to it...")
      dims <- seq(min(dims), ncol(reduction), 1)
    }
    reduction <- reduction[, dims]
  }
  
  umap.embeddings <- uwot::umap(X = reduction,
                                n_neighbors = as.integer(x = n.neighbors),
                                n_components = as.integer(x = n.components),
                                metric = metric,
                                n_epochs = n.epochs,
                                learning_rate = learning.rate,
                                min_dist = min.dist,
                                spread = spread,
                                set_op_mix_ratio = set.op.mix.ratio,
                                local_connectivity = local.connectivity,
                                repulsion_strength = repulsion.strength,
                                negative_sample_rate = negative.sample.rate,
                                a = a, b = b,
                                fast_sgd = uwot.sgd,
                                verbose = verbose)
  rownames(umap.embeddings) <- rownames(reduction)
  return(umap.embeddings)
}


plotCompositionBarplots <- function(groups, sample.factor=NULL, 
                                    show.entropy=TRUE,show.size=TRUE, show.composition=TRUE,legend.height=0.2) {
  ## param checking
  if (is.null(groups)) {
    stop('groups factor on the cells needs to be specified')
  }
  
  groups <- as.factor(groups)
 
  xt <- table(sample.factor[match(names(groups),names(sample.factor))],groups)
  xt <- xt[rowSums(xt)>0,]; xt <- xt[,colSums(xt)>0]
  
  df <- reshape2::melt(xt); colnames(df) <- c("sample","cluster","f");  df$f <- df$f/colSums(xt)[as.character(df$cluster)]
  clp <- ggplot2::ggplot(df, ggplot2::aes(x=factor(cluster, levels=levels(groups)),y=f,fill=sample)) +
    ggplot2::geom_bar(stat='identity') + ggplot2::xlab('cluster') + ggplot2::ylab('fraction of cells') + ggplot2::theme_bw() +
    ggplot2::scale_y_continuous(expand=c(0, 0))
  
  if(!show.size && !show.entropy)
    return(clp);
  
  # extract legend
  leg <- cowplot::get_legend(clp + ggplot2::theme(legend.position="bottom"))
  pl <- list(clp + ggplot2::theme(legend.position="none"));
  
  if(show.entropy) {
    if (!requireNamespace("entropy", quietly=T))
      stop("You need to install 'entropy' package to use 'show.entropy=T'")
    
    n.samples <- nrow(xt);
    ne <- 1-apply(xt, 2, entropy::KL.empirical, y2=rowSums(xt), unit=c('log2')) / log2(n.samples) # relative entropy
    enp <- ggplot2::ggplot(data.frame(cluster=factor(colnames(xt),levels=levels(groups)),entropy=ne), ggplot2::aes(cluster, entropy)) +
      ggplot2::geom_bar(stat='identity',fill='grey65') + ggplot2::ylim(0,1) +
      ggplot2::geom_hline(yintercept=1, linetype="dashed", color = "grey30") + ggplot2::theme_bw()
    pl <- c(pl,list(enp))
  }
  
  if(show.size) {
    szp <- ggplot2::ggplot(data.frame(cluster=factor(colnames(xt),levels=levels(groups)), cells=colSums(xt)), ggplot2::aes(cluster,cells)) +
      ggplot2::geom_bar(stat='identity') + ggplot2::scale_y_continuous(trans='log10') + ggplot2::theme_bw() + ggplot2::ylab('number of cells')
    pl <- c(pl,list(szp))
  }
  
  return(pl)
}


getCorrected <- function(y, mod, svs) {
  X = cbind(mod, svs)
  Hat = solve(t(X) %*% X) %*% t(X)
  beta = (Hat %*% t(y))
  P = ncol(mod)
  return(y - t(as.matrix(X[,-c(1:P)]) %*% beta[-c(1:P),]))
}

#' differential expression based on DESeq
#' @mat speudo-bulk count matrix - rows are genes columns are sample
#' @meta sample meta-table, the first column contains condition need to be compare
#' @batch.column which column contains batch
DEseq <- function(mat, meta=NULL, batch.column=NULL, pca.plot=TRUE, test.name="Wald",
                  target.level="AKI", ref.level="Ref"){
  if(is.null(meta)){
    stop("please provide meta-table contains sample information")
  }
  design.formula <- as.formula(paste('~ ', paste(colnames(meta)[c(1)], collapse=' + ')))
  
  if(is.null(batch.column)){
    cnts.norm <- mat  %>% DESeq2::DESeqDataSetFromMatrix(colData = meta, design= design.formula)  %>%
      DESeq2::estimateSizeFactors()  %>% DESeq2::counts(normalized=TRUE)
  } else{
    design.formula.batch <- as.formula(paste('~ ', paste(colnames(meta)[c(1, batch.column)], collapse=' + ')))
    #design.formula.2 <- as.formula(paste('~ ', paste(colnames(meta)[c(batch.column)], collapse=' + ')))
    #dds <- DESeqDataSetFromMatrix(countData=mat, colData=meta, design = design.formula) %>% DESeq()
    
    # original norm matrix
    cnts.norm <- mat %>% DESeq2::DESeqDataSetFromMatrix(colData = meta, design= design.formula.batch)  %>%
      DESeq2::estimateSizeFactors()  %>% DESeq2::counts(normalized=TRUE)
    
    # remove batch in the normalized gene matrix
    mod <- model.matrix(design.formula, data=meta)
    counts_batch_corrected <- getCorrected(cnts.norm, mod, meta[, batch.column])
    design.formula <- design.formula.batch
  }
  
  if(is.null(batch.column)){
    cnts <- cnts.norm
    counts_batch_corrected <- NULL
  } else{
    cnts <- counts_batch_corrected
  }
  
  # remove gene with 0 expression
  cnts <- cnts[rowSums(cnts) > 0, ]
  cnts <- cnts[, colSums(cnts) >0]
  
  if(pca.plot){
    cnts.pca <- prcomp(t(cnts), center = TRUE, scale. =TRUE)
    sample.pc <- as.data.frame(cnts.pca$x[, 1:2])
    sample.pc$condition <- meta[, 1]
    sample.pc$batch <- meta[, 2]
    
    ## plot
    p1 <- ggplot(sample.pc, aes(x=PC1, y=PC2, color=condition)) + geom_point() + theme_bw()
    p2 <- ggplot(sample.pc, aes(x=PC1, y=PC2, color=batch)) + geom_point() + theme_bw()
    p3 <- ggplot(sample.pc, aes(x=PC1, y=batch)) + geom_point(size=2, alpha=0.6) + xlab("PC1") +
      ylab("batch") + stat_cor(method = "pearson") + theme_bw()
    print(cowplot::plot_grid(plotlist=list(p1, p2, p3), ncol=2))
  }
  
  
  # differential expression test
  if(test.name == "Wald"){
    res1 <- DESeq2::DESeqDataSetFromMatrix(mat, meta, design=design.formula) %>%
      DESeq2::DESeq(quiet=T, test=test.name) %>%
      DESeq2::results(contrast=c(colnames(meta)[1], target.level, ref.level),
                      cooksCutoff = FALSE,
                      independentFiltering = TRUE) %>% as.data.frame
    res1 <- res1[order(res1$padj, decreasing=F), ]
  } 
  
  return(list(norm=cnts.norm, corrected=counts_batch_corrected, diffgene=res1, model=design.formula))
}

#' fisher-exact statistics for enrichment of dev markers in clusters
#' @param gene.set - dev markers
fTest <- function(sample.de.genes, all.de.genes, gene.set){
  n <- length(intersect(sample.de.genes, gene.set))
  nc.genes <- all.de.genes[-which(all.de.genes %in% sample.de.genes)]
  no <- length(intersect(nc.genes, gene.set))
  
  c.table <- matrix(c(n, length(sample.de.genes)-n, no, length(all.de.genes)-no), nrow=2)
  est <- fisher.test(c.table)
  #c.table <- matrix(c(n.all.genes - length(union(de.genes, gene.set)), length(de.genes) - n.de.in.gs,
  #                    length(gene.set) - n.de.in.gs, n.de.in.gs), nrow = 2)
  stat <- data.frame(pval=est$p.value, odds.ratio=est$estimate, conf_low=est$conf.int[1], conf_high=est$conf.int[2])
  return(stat)
}

#' annotate peak regions in human
#' @param vector of peak lists (e.g. chr1-180738-180958)
#' @return data.frame contains annotated features for each region
annotPeaks <- function(peak.list){
  require(ChIPseeker)
  require(org.Hs.eg.db)
  require(TxDb.Hsapiens.UCSC.hg38.knownGene)
  
  peak.df <- data.frame(chr=peak.list) %>%
    tidyr::separate(col=chr, into = c("chr", "start", "end"), sep = "-")
  peak.gr <- with(peak.df, GRanges(as.character(chr),
                                   IRanges(as.numeric(as.character(start)), as.numeric(as.character(end)))))
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  peak_anno <- annotatePeak(peak.gr, tssRegion = c(-5000, 5000),
                            TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene, annoDb ="org.Hs.eg.db")
  
  peaksInfo <- as.data.frame(peak_anno@anno)
  regionsAnnot <- as.data.frame(peak_anno@detailGenomicAnnotation)
  # transform to numeric value
  regionsAnnot <- sapply(regionsAnnot, as.numeric)
  # add to peaksInfo
  peaksInfo <- cbind(peaksInfo, regionsAnnot)
  return(peaksInfo)
}

fTest <- function(sample.de.genes, all.de.genes, gene.set){
  n <- length(intersect(sample.de.genes, gene.set))
  nc.genes <- all.de.genes[-which(all.de.genes %in% sample.de.genes)]
  no <- length(intersect(nc.genes, gene.set))
  
  c.table <- matrix(c(n, length(sample.de.genes)-n, no, length(all.de.genes)-no), nrow=2)
  est <- fisher.test(c.table)
  #c.table <- matrix(c(n.all.genes - length(union(de.genes, gene.set)), length(de.genes) - n.de.in.gs,
  #                    length(gene.set) - n.de.in.gs, n.de.in.gs), nrow = 2)
  stat <- data.frame(pval=est$p.value, odds.ratio=est$estimate, conf_low=est$conf.int[1], conf_high=est$conf.int[2])
  return(stat)
}


#' Compute relationship between arbitrary regions and genes using either Hi-C or coaccessibility
#'
#' @param regions Genomic regions in granges format
#' @param conns Dataframe of peak to peak Hi-C/coaccessibility
#' @param link.promoter Include peaks in gene promoters
#' @param promoter.region Specify the window around the TSS that counts as a promoter
#' @param anno.level Specify "gene" or "transcript" for a gene/transcript level annotation
#' @param region.name Column name of region identifier. Defaults to peak names
#' @param weight.col Column name of weights (i.e. when using coaccessibility). Defaults to 1
#'
#' @return Matrix of region to gene Hi-C/coaccessibility connections
#'
#' @import GenomicRanges
#' @import ChIPseeker
#' @import Matrix
#' @export
#'
RegionGeneLinks <- function(regions,
                            conns,
                            link.promoter = T,
                            promoter.region = c(-3000, 3000),
                            anno.level = "transcript",
                            region.name = NULL,
                            weight.col = NULL){
  require(org.Hs.eg.db)
  require(TxDb.Hsapiens.UCSC.hg38.knownGene)
  
  ## Get region names
  regions <- unique(regions)
  regions.peaks <- granges2peak(regions)
  if (is.null(region.name)) {
    names(regions) <- regions.peaks
  } else {
    names(regions) <- regions@elementMetadata[[region.name]]
  }
  
  ## Annotate regions
  regions.anno <- annotatePeak(regions, tssRegion = promoter.region, level = anno.level,
                               TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,
                               annoDb = "org.Hs.eg.db")
  regions.anno <- as.data.frame(regions.anno)
  
  ## Find overlaps between peak1 and regions
  peak1.gr <- peak2granges(conns$Peak1)
  regions.conns.ix <- findOverlaps(peak1.gr, regions)
  
  ## Annotate peak2
  peak2.gr <- peak2granges(conns$Peak2)
  peak2.gr.anno <- annotatePeak(peak2.gr, tssRegion = promoter.region, level = anno.level,
                                TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,
                                annoDb = "org.Hs.eg.db")
  peak2.gr.anno <- peak2.gr.anno@anno
  
  ## Get weights (otherwise set to 1 if null)
  if (!is.null(weight.col)) {
    stopifnot(weight.col %in% colnames(conns))
    peak2.gr.anno$weights <- conns[[weight.col]]
  } else {
    peak2.gr.anno$weights <- 1
  }
  
  regions.conns.ix.vector <- as.character(regions.conns.ix@to)
  names(regions.conns.ix.vector) <- as.character(regions.conns.ix@from)
  regions.conns.ix.list <- lapply(unique(regions.conns.ix.vector), function(x) {
    names(regions.conns.ix.vector[regions.conns.ix.vector == x])
  })
  names(regions.conns.ix.list) <- names(regions)[as.integer(unique(regions.conns.ix.vector))]
  
  ## Link region to a gene if the region is in the promoter region
  if (link.promoter) {
    is.promoter <- grepl("Promoter", regions.anno$annotation)
    region.gene.weights.list <- lapply(1:length(regions), function(i) {
      if (is.promoter[[i]]){
        w <- 1
        names(w) <- regions.anno[["SYMBOL"]][[i]]
      } else {
        w <- c()
      }
      w
    })
    names(region.gene.weights.list) <- names(regions)
  } else {
    region.gene.weights.list <- lapply(1:length(regions), function(i) { c() })
    names(region.gene.weights.list) <- names(regions)
  }
  
  ## Link coaccessibility
  for (h in names(regions.conns.ix.list)) {
    peaks.ix <- as.integer(regions.conns.ix.list[[h]])
    peaks.anno.gr <- peak2.gr.anno[peaks.ix]
    peaks.anno.gr <- peaks.anno.gr[grepl("Promoter", peaks.anno.gr$annotation)]
    if (length(peaks.anno.gr) > 0) {
      w <- peaks.anno.gr$weights
      names(w) <- peaks.anno.gr$SYMBOL
      region.gene.weights.list[[h]] <- c(region.gene.weights.list[[h]], w)
    }
  }
  
  region.gene.weights.list <- lapply(region.gene.weights.list, function(w) {
    if (length(w) > 1) w <- tapply(w, names(w), max)
    w[!(is.na(w) | is.na(names(w)))]
  })
  region.gene.weights.list <- region.gene.weights.list[sapply(region.gene.weights.list, function(x) !is.null(x))]
  region.gene.weights.list <- region.gene.weights.list[sapply(region.gene.weights.list, function(x) length(x) > 0)]
  region.gene.weights.df <- lapply(names(region.gene.weights.list), function(i) {
    weights <- region.gene.weights.list[[i]]
    data.frame(region = i, gene = names(weights), weight = weights,
               stringsAsFactors = F)
  })
  region.gene.weights.df <- do.call(rbind, region.gene.weights.df)
  rownames(region.gene.weights.df) <- NULL
  return(region.gene.weights.df)
}

#' Function for converting a granges object into a character vector of peak names
#'
#' @param gr GenomicRanges object
#' @param delim Peak name delimiter
#'
#' @return Character vector of peak names
#' @import GenomicRanges
#' @export
#'
granges2peak <- function(gr, delim = c(":", "-")) {
  paste0(seqnames(gr), delim[[1]],
         start(gr), delim[[2]],
         end(gr))
}

### sample-level clustering functions
fitLogistic <- function(exp, cond, comparison=c("AKI", "CKD")){
  # exp: expression value
  # cond: condition
  cond.index <- which(cond %in% comparison)
  exp <- exp[cond.index]
  cond <- cond[cond.index]
  cond <- as.factor(cond)
  
  model <- glm((cond)~exp, family = binomial(link = 'logit'))
  
  ## caculate AUC
  log_predict <- predict(model, newdata=as.data.frame(exp), type = "response")
  log_predict <- ifelse(log_predict > 0.5,1,0)
  
  # binarize
  cond.bin <- rep(0, length(cond))
  cond.bin[cond==comparison[1]] <- 0
  cond.bin[cond==comparison[2]] <- 1
  
  auc <- pROC::auc(cond.bin, log_predict)
  
  return(list(model, auc))
}

getModels <- function(patient_rpkm_genesets, patients){
  ### train classifier for each gene
  
  aucs <- list()
  models <- list()
  for(i in rownames(patient_rpkm_genesets)){
    exp <- patient_rpkm_genesets[i, ]
    cond <- patients$condition
    res <- fitLogistic(exp, cond, comparison=c("AKI", "CKD"))
    aucs[[i]] <- data.frame(gene=i, auc=res[[2]])
    models[[i]] <- res[[1]]
  }
  names(models) <- rownames(patient_rpkm_genesets)
  
  aucs <- do.call(rbind, aucs)
  aucs <- aucs[order(aucs$auc, decreasing = T), ]
  return(list(aucs, models))
}

#' plot fitted value from logistic regression
#' @param model gene specific model object
plotLogitFit <- function(gene.model, colorvec=NULL){
  fitted.values <- gene.model[[1]]$fitted.values
  if(!is.null(colorvec)){
    plot(sort(fitted.values), col=colorvec[order(fitted.values)], pch=20, xlab="", ylab="", ylim=c(0, 1))
    abline(h=0.5, lty=2, col="grey")
    text(x=3, y=0.93, names(gene.model))
  } else{
    plot(sort(fitted.values), pch=20, xlab="", ylab="", ylim=c(0, 1))
    abline(h=0.5, lty=2, col="grey")
  }
}

#' plot fitted value from logistic regression
#' @param model gene specific model object
plotLogitFit2inte <- function(gene.model1, gene.model2, colorvec=NULL){
  fitted.values1 <- gene.model1[[1]]$fitted.values
  fitted.values2 <- gene.model2[[1]]$fitted.values
  fitted.values <- c(fitted.values1, fitted.values2)
  if(!is.null(colorvec)){
    plot(sort(fitted.values), col=colorvec[order(fitted.values)], pch=20, xlab="", ylab="", ylim=c(0, 1))
    abline(h=0.5, lty=2, col="grey")
    text(x=3, y=0.93, names(gene.model1))
  } else{
    plot(sort(fitted.values), pch=20, xlab="", ylab="", ylim=c(0, 1))
    abline(h=0.5, lty=2, col="grey")
  }
}

# project stratified gene expression to tsne plot
plotTsneGene <- function(exp, tSNE, cond.comp, out=1, gene.name=NULL){
  
  cond.comp <- as.factor(cond.comp)
  
  model <- glm((cond.comp)~tSNE$Y,family="binomial")
  slope <- coef(model)[2]/(-coef(model)[3])
  intercept <- coef(model)[1]/(-coef(model)[3]) 
  
  # scale
  exp[exp>quantile(exp, 0.7)] <- quantile(exp, 0.7)
  exp[exp<quantile(exp, 0.3)] <- quantile(exp, 0.3)
  
  grad <- (exp - min(exp)) / (max(exp) - min(exp))
  plot(tSNE$Y[,1],tSNE$Y[,2],cex=1.5,las=1,xaxt="n",yaxt="n",
       col=colorRampPalette(c("grey90","grey80","navy"))(50)[grad*49+1],
       pch=20,line=.5,xlab="tSNE1", ylab="tSNE2",xlim=c(min(tSNE$Y[,1])*out,max(tSNE$Y[,1])*out),
       ylim=c(min(tSNE$Y[,2])*out,max(tSNE$Y[,2])*out), main=gene.name)
  axis(side = 1, labels = FALSE, tck = -0.01);axis(side = 2, labels = FALSE, tck = -0.01)
  abline(intercept,slope,lty=2)
}

plotContourTsne <- function(tSNE, cond.comp, title=NULL, col, cex=2.5){
  cond.comp <- as.factor(cond.comp)
  
  model <- glm((cond.comp)~tSNE$Y,family="binomial")
  slope <- coef(model)[2]/(-coef(model)[3])
  intercept <- coef(model)[1]/(-coef(model)[3]) 
  
  
  z <- kde2d(tSNE$Y[,1], tSNE$Y[,2], n=100, lims = c(min(tSNE$Y[,1])*out, max(tSNE$Y[,1])*out, min(tSNE$Y[,2])*out, max(tSNE$Y[,2])*out))
  plot(tSNE$Y[,1],tSNE$Y[,2],type="n",xaxt="n",yaxt="n",line=.5,xlab="tSNE1", ylab="tSNE2",xlim=c(min(tSNE$Y[,1])*out,max(tSNE$Y[,1])*out),ylim=c(min(tSNE$Y[,2])*out,max(tSNE$Y[,2])*out),main=title)
  #contour(z, drawlabels=FALSE, nlevels=5, col=c("red","grey","blue2","orange3","firebrick"), add=TRUE)
  contour(z, drawlabels=FALSE, nlevels=5, col=c("grey"), add=TRUE)
  points(tSNE$Y[,1], tSNE$Y[,2],cex=cex, pch=20, col=col)
  abline(intercept,slope,lty=2, col="lightgrey")
}

plotOverlapSets <- function(genesets1, genesets2){
  genesets.union <- unique(c(genesets1, genesets2))
  sc.vec <- rep(0, length(genesets.union))
  sn.vec <- rep(0, length(genesets.union))
  
  df <- data.frame(sn=sn.vec, sc=sc.vec)
  df$genesets <- genesets.union
  df$sn[df$genesets %in% genesets1] <- 1
  df$sc[df$genesets %in% genesets2] <- 1
  
  d <- vennCounts(df[, 1:2])
  vennDiagram(d, circle.col=c("#507EB3", "#4D7880"))
}
