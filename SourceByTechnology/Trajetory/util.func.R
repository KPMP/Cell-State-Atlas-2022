
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



