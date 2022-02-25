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

#' plot significant heatmap from LR analysis - customized from CellChat
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

