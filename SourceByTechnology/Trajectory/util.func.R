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


#' plot gene expression along with psedo-time
#' 
plot_gene_psedotime <- function(mat, gene.name=NULL, pseudotime, groups, dot.size=1, alpha=0.05,
                                clpalatte=NULL, difftest.plot=TRUE, spline.df=3,
                                condition.smooth=FALSE){
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
  
  if(condition.smooth){
    q <- ggplot(aes(pseudotime, exp, color=annot), data = df) + 
      geom_point(size=dot.size, alpha=alpha, data=df) + 
      geom_smooth(se = FALSE, method = "gam", formula = y ~ s(log(x))) + theme_bw() +
      scale_color_manual(values=clpalatte)
  }
  
  return(q)
}


rescale.and.center <- function(x, center=F, max.quantile=0.99) {
  mx <- quantile(abs(x),max.quantile) # absolute maximum
  if(mx==0) mx<-max(abs(x)) # in case the quantile squashes all the signal
  x[x>mx] <- mx; x[x< -1*mx] <- -1*mx; # trim
  if(center) x <- x-mean(x) # center
  x/max(abs(x)); # scale
}

#' plot heatmap of smoothed gene expression along psedotime
#' @param gene_mat smoothed gene expression matrix, rows are genes columns are cells
plot_trajectory_heatmap <- function(gene_mat, cellannot, col.row=NULL, col.df=NULL, colpal, row.font.size=8){
  
  #pheatmap::pheatmap(gene_mat, cluster_rows=F, cluster_cols=F, show_colnames=F,
  #                   annotation_row=row.df, annotation_col=col.df, fontsize_row=1)
  if(!is.null(colpal)){
    col.pt <- colpal[names(colpal) %in% names(cellannot)]
    col.pt <- col.pt[match(rownames(col.df), names(col.pt))]
    names(col.pt) <- col.df[, 1]
  }
  
  if(!is.null(col.df)){
    ha <- ComplexHeatmap::HeatmapAnnotation(celltype=col.df[, 1],
                                            which='column', show_annotation_name=FALSE, show_legend=T, border=T, 
                                            col=list(celltype=col.pt))
  }
  
  if(!is.null(col.row)){
    hl <- ComplexHeatmap::HeatmapAnnotation(cluster=names(col.row),
                                            which='row', show_annotation_name=FALSE, show_legend=T, border=T, 
                                            col=list(cluster=col.row))
  }
  
  if(is.null(col.df)){
    hm <- ComplexHeatmap::Heatmap(t(scale(t(gene_mat))), cluster_columns = F, cluster_rows = F,show_column_names = F,border=T, 
                                  name='expression', show_heatmap_legend = F, show_row_dend = F, 
                                  show_row_names=T, left_annotation=hl,  row_names_gp = gpar(fontsize = row.font.size))
  } else if(is.null(col.row)){
    hm <- ComplexHeatmap::Heatmap(t(scale(t(gene_mat))), cluster_columns = F, cluster_rows = F,show_column_names = F,border=T, 
                                  name='expression', show_heatmap_legend = F, show_row_dend = F, 
                                  show_row_names=T, top_annotation = ha,  row_names_gp = gpar(fontsize = row.font.size))
  } else{
    hm <- ComplexHeatmap::Heatmap(t(scale(t(gene_mat))), cluster_columns = F, cluster_rows = F,show_column_names = F,border=T, 
                                  name='expression', show_heatmap_legend = F, show_row_dend = F, 
                                  show_row_names=T, top_annotation = ha, left_annotation=hl,
                                  row_names_gp = gpar(fontsize = row.font.size))
  }
  
  #labeled.genes <- rownames(gene_mat)
  #labeled.genes <- c(labeled.genes, c("AQP1", "DCDC2", "HAVCR1", "LCN2", "PLSCR1", "PROM1", "UMOD", "VCAM1", "ITGB6", "EGF", "SLC12A1"))
  #hm + ComplexHeatmap::rowAnnotation(link = ComplexHeatmap::anno_mark(at = match(labeled.genes, rownames(gene_mat)), 
  #                                                                    labels = labeled.genes))
  hm
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

#' plot summary of DE genes
#' @param cao list of cacoa object
plotDEsummary <- function(cao, p.adjust = TRUE, pvalue.cutoff = 0.05, 
                           show.jitter = FALSE, jitter.alpha = 0.05, palette=cao[[1]]$cell.groups.palette,
                           type = "bar", notch = TRUE, show.whiskers = TRUE, show.regression = TRUE){
  
  # turn resampling option off
  show.resampling.results = FALSE
  
  df <- list()
  for(i in 1:length(cao)){
    de.raw <- cao[[i]]$test.results$deFixed
    
    if (show.resampling.results) {
      if (!all(unlist(lapply(de.raw, function(x) !is.null(x$subsamples))))) {
        warning("resampling results are missing for at least some cell types, falling back to point estimates. Please rerun estimatePerCellTypeDE() with resampling='bootstrap' or resampling='loo'")
        rl <- lapply(de.raw, `[[`, "res")
      }
      else {
        subsamples <- lapply(de.raw, `[[`, "subsamples")
        rl <- unlist(subsamples, recursive = FALSE) %>% setNames(rep(names(de.raw), 
                                                                     sapply(subsamples, length)))
      }
    }
    else {
      rl <- lapply(de.raw, `[[`, "res")
    }
    
    df[[i]] <- do.call(rbind, lapply(1:length(rl), function(j){
      if (p.adjust) {
        ndiff <- sum(na.omit(rl[[j]]$padj <= pvalue.cutoff))
      }
      else {
        ndiff <- sum(na.omit(rl[[j]]$pvalue <= pvalue.cutoff))
      }
      data.frame(Type = names(rl)[j], value = ndiff, stringsAsFactors = FALSE)
    }))
    
  }
  df <- df %>% Reduce(function(dtf1 ,dtf2) left_join(dtf1, dtf2,by="Type"), .)
  df$value <- rowSums(df[, 2:ncol(df)])
  

  p <- cacoa:::plotCellTypeSizeDep(df, cao[[1]]$cell.groups, palette = palette, 
                             ylab = "number of DE genes", yline = NA, show.whiskers = show.whiskers, 
                             show.regression = show.regression, plot.theme=cao[[1]]$plot.theme)

  return(p)
}


#' plot Go Heatmap from cao object
plotOntologyHeat <- function (cao, genes = "up", type = "GO", subtype = "BP", min.genes = 1, 
                              p.adj = 0.05, legend.position = "left", selection = "all", 
                              n = 20, clusters = TRUE, cluster.name = NULL, cell.subgroups = NULL, 
                              color.range = NULL, palette = NULL, row.order = TRUE, col.order = TRUE, 
                              legend.title = NULL, row.dendrogram = FALSE, col.dendrogram = FALSE, up.qantile=0.75){
  if (!is.null(cell.subgroups) && (length(cell.subgroups) == 1)) 
    stop("'cell.subgroups' must contain at least two groups. Please use plotOntology instead.")
  if (is.null(selection) || (!selection %in% c("unique", "common", "all"))) 
    stop("'selection' must be one of the following: 'unique', 'common', or 'all'.")
  
  name <- if (is.null(cluster.name)) 
      cacoa:::getOntClustField(subtype, genes)
    else cluster.name
    
  if (is.null(cao[[1]]$test.results[[type]][[name]])) {
      if (!is.null(cluster.name)) 
        stop("Can't find the results for ", cluster.name)
  }
  
  ont.sum.list <- lapply(cao, function(r){
    ont.sum <- r$test.results[[type]][[name]]$df %>% cacoa:::groupOntologiesByCluster(field = "ClusterName") %>% 
      {
        . * -1
      } %>% as.data.frame()
    ont.sum$term <- rownames(ont.sum)
    ont.sum
  })
  
  ## merge list
  ont.sum <- merge(ont.sum.list[[1]], ont.sum.list[[2]], by=c("term"), all=TRUE)
  ont.sum[is.na(ont.sum)] <- 0
  colnames(ont.sum) <- gsub("[.]x|[.]y", "", colnames(ont.sum))
  terms <- ont.sum$term
  
  ont.sum <- sapply(split.default(ont.sum[, -1], names(ont.sum)[-1]), rowMeans, na.rm = TRUE)
  rownames(ont.sum) <- terms
  
  if (nrow(ont.sum) == 0) {
    warning(paste0("No ontologies pass the filtration for type=", 
                   type, ", subtype=", subtype, " and genes=", genes))
    return(ggplot())
  }
  if (selection == "unique") {
    ont.sum %<>% .[rowSums(abs(.) > 0) == 1, , drop = FALSE]
  }
  else if (selection == "common") {
    ont.sum %<>% .[rowSums(abs(.) > 0) > 1, , drop = FALSE]
  }
  if (nrow(ont.sum) == 0) 
    stop("Nothing to plot. Try another selection.")
  tmp <- as.matrix(ont.sum %>% .[order(rowSums(.), decreasing = TRUE), 
  ] %>% .[1:pmin(nrow(.), n), ]) %>% {
    .[, !colSums(.) == 0]
  }
  
  tmp %<>% .[, order(apply(tmp, 2, function(r){ sum(r > 0 ) }), decreasing=T)]
  ## scale
  up.value <- quantile(tmp, up.qantile)
  tmp[tmp > up.value] <- up.value
  
  if (is.null(color.range)) {
    color.range <- c(min(0, min(tmp, na.rm = TRUE)), max(tmp, 
                                                         na.rm = TRUE))
    if (color.range[2] > 20) {
      warning("Shrinking minimum adj. P value to -log10(20) for plotting.")
      color.range[2] <- 20
    }
    tmp %<>% pmax(color.range[1]) %>% pmin(color.range[2])
    title = "-log10(adj. P)"
  }
  else {
    if (is.null(legend.title)) 
      title = "Bin"
    else title = legend.title
  }
  pal <- if (genes == "up") {
    circlize::colorRamp2(c(color.range[1], color.range[2]), 
                         c("grey98", "red"))
  }
  else if (genes == "down") {
    circlize::colorRamp2(c(color.range[1], color.range[2]), 
                         c("grey98", "blue"))
  }
  else {
    circlize::colorRamp2(c(color.range[1], color.range[2]), 
                         c("grey98", "darkgreen"))
  }
  ComplexHeatmap::Heatmap(tmp, col = pal, border = TRUE, show_row_dend = row.dendrogram, 
                          cluster_columns=FALSE, heatmap_legend_param = list(title = title), 
                          row_names_max_width = unit(8, "cm"), row_names_gp = grid::gpar(fontsize = 10))
}

#' plot cell-cell communication probability
#' @param color.use factor specifies colors for each celltype
plot_communication <- function(obj, module, cutoff=0.1, rescale=FALSE, feature='weight', norm=T, color.use=NULL,
                               cell.order=NULL){
  if(feature == 'weight'){
    mat <- obj@net$weight
  } else if(feature=='count'){
    mat <- obj@net$count
    if(norm){
      mat <- t(apply(mat, 1, function(r) r/sum(r)))
    }
  }
  
  
  #par(mfrow = c(3,4), xpd=TRUE)
  row.idx <- which(rownames(mat) %in% module)
  weight.matrix <- list()
  
  for (i in 1:length(row.idx)) {
    # remove other modules
    mat2 <- matrix(0, nrow = nrow(mat)-length(row.idx)+1, ncol = ncol(mat)-length(row.idx)+1)
    rownames(mat2) <- rownames(mat)[-row.idx[-i]]
    colnames(mat2) <- colnames(mat)[-row.idx[-i]]
    
    mat2[rownames(mat2) %in% rownames(mat)[row.idx[i]], ] <- mat[row.idx[i], -row.idx[-i]]
    
    
    if(rescale){
      values <- mat[row.idx[i], -row.idx[-i]]
      mat2[rownames(mat2) %in% rownames(mat)[row.idx[i]], ] <- values/(max(values) - min(values))
    }
    
    mat2[mat2<cutoff] <- 0
    
    groupSize <- as.numeric(table(obj@idents))[-row.idx[-i]]
    
    # remove self interaction
    diag(mat2)  <- 0
    
    #if(use.module.color){
    #  color.use[which(rownames(mat2) %in% rownames(mat)[row.idx[i]])] <- rownames(mat)[row.idx[i]]
    #}
    
    ## re-order
    if(!is.null(cell.order)){
      order <- rownames(mat2)[match(cell.order, rownames(mat2))] 
      order <- order[!is.na(order)]
      mat2 <- mat2[order, order]
      
    } else{
      order <- rownames(mat2)
    }
    
    if(is.null(color.use)){
      color = scPalette(nrow(mat2))
    } else{
      color <- color.use[match(order, names(color.use))]
    }
    
    
    netVisual_circle(mat2, vertex.weight = groupSize, color.use=color, 
                     weight.scale = T, edge.weight.max = max(mat2), 
                     title.name = rownames(mat)[row.idx[i]])
    
    weight.matrix[[i]] <- mat2
  }
  names(weight.matrix) <- module
  return(weight.matrix)
}

#' plot heatmap of significant LR pairs
plotSigHeatmap <- function (object, signaling = NULL, pattern = c("outgoing", "incoming", "all"), 
                            slot.name = "netP", color.use = NULL, color.heatmap = "BuGn", 
                            title = NULL, width = 10, height = 8, font.size = 8, font.size.title = 10, 
                            cluster.rows = FALSE, cluster.cols = FALSE, sel.group=NULL){
  
  pattern <- match.arg(pattern)
  if (length(slot(object, slot.name)$centr) == 0) {
    stop("Please run `netAnalysis_computeCentrality` to compute the network centrality scores! ")
  }
  centr <- slot(object, slot.name)$centr
  outgoing <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
  incoming <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
  dimnames(outgoing) <- list(levels(object@idents), names(centr))
  dimnames(incoming) <- dimnames(outgoing)
  for (i in 1:length(centr)) {
    outgoing[, i] <- centr[[i]]$outdeg
    incoming[, i] <- centr[[i]]$indeg
  }
  if(pattern == "outgoing"){
    mat <- t(outgoing)
    legend.name <- "Outgoing"
  } else if (pattern == "incoming") {
    mat <- t(incoming)
    legend.name <- "Incoming"
  } else if (pattern == "all"){
    mat <- t(outgoing + incoming)
    legend.name <- "Overall"
  }
  
  if(!is.null(sel.group)){
    mat <- mat[, sel.group]
  }
  
  # remove 0 rows
  if(length(which(rowSums(mat)==0))>0){
    mat <- mat[-which(rowSums(mat)==0), ]
  }
  
  if (is.null(title)) {
    title <- paste0(legend.name, " signaling patterns")
  } else {
    title <- paste0(paste0(legend.name, " signaling patterns"), 
                    " - ", title)
  }
  
  if (!is.null(signaling)) {
    mat1 <- mat[rownames(mat) %in% signaling, , drop = FALSE]
    mat <- matrix(0, nrow = length(signaling), ncol = ncol(mat))
    idx <- match(rownames(mat1), signaling)
    mat[idx[!is.na(idx)], ] <- mat1
    dimnames(mat) <- list(signaling, colnames(mat1))
  }
  
  mat.ori <- mat
  mat <- sweep(mat, 1L, apply(mat, 1, max), "/", check.margin = FALSE)
  
  mat[mat == 0] <- NA
  
  if (is.null(color.use)) {
    color.use <- scPalette(length(colnames(mat)))
  }
  color.heatmap.use = (grDevices::colorRampPalette((RColorBrewer::brewer.pal(n = 9, 
                                                                             name = color.heatmap))))(100)
  
  df <- data.frame(group = colnames(mat))
  rownames(df) <- colnames(mat)
  names(color.use) <- colnames(mat)
  col_annotation <- HeatmapAnnotation(df = df, col = list(group = color.use), 
                                      which = "column", show_legend = FALSE, show_annotation_name = FALSE, 
                                      simple_anno_size = grid::unit(0.2, "cm"))
  
  ha2 = ComplexHeatmap::HeatmapAnnotation(Strength = anno_barplot(colSums(mat.ori), 
                                                                  border = FALSE, gp = gpar(fill = color.use, col = color.use)), 
                                          show_annotation_name = FALSE)
  
  pSum <- rowSums(mat.ori)
  pSum.original <- pSum
  pSum <- -1/log(pSum)
  pSum[is.na(pSum)] <- 0
  idx1 <- which(is.infinite(pSum) | pSum < 0)
  if (length(idx1) > 0) {
    values.assign <- seq(max(pSum) * 1.1, max(pSum) * 1.5, 
                         length.out = length(idx1))
    position <- sort(pSum.original[idx1], index.return = TRUE)$ix
    pSum[idx1] <- values.assign[match(1:length(idx1), position)]
  }
  
  ha1 = rowAnnotation(Strength = anno_barplot(pSum, border = FALSE), 
                      show_annotation_name = FALSE)
  
  if (min(mat, na.rm = T) == max(mat, na.rm = T)) {
    legend.break <- max(mat, na.rm = T)
  } else {
    legend.break <- c(round(min(mat, na.rm = T), digits = 1), 
                      round(max(mat, na.rm = T), digits = 1))
  }
  
  ht1 = Heatmap(mat, col = color.heatmap.use, na_col = "white", 
                name = "Relative strength", bottom_annotation = col_annotation, 
                top_annotation = ha2, right_annotation = ha1, cluster_rows = cluster.rows, 
                cluster_columns = cluster.rows, row_names_side = "left", 
                row_names_rot = 0, row_names_gp = gpar(fontsize = font.size), 
                column_names_gp = gpar(fontsize = font.size), width = unit(width, 
                                                                           "cm"), height = unit(height, "cm"), column_title = title, 
                column_title_gp = gpar(fontsize = font.size.title), column_names_rot = 90, 
                heatmap_legend_param = list(title_gp = gpar(fontsize = 8, 
                                                            fontface = "plain"), title_position = "leftcenter-rot", 
                                            border = NA, at = legend.break, legend_height = unit(20, 
                                                                                                 "mm"), labels_gp = gpar(fontsize = 8), grid_width = unit(2, 
                                                                                                                                                          "mm")))
  return(list(heatmap=ht1, mat=mat))
}
