require(Giotto)
require(ggplot2)
require(Matrix)
require(patchwork)

# A function for adding RCTD weights without removing beads 
readPuck <- function(puck, puck.addr, allweights, instrs){
  ssg.puck <- readSlideSeqGiotto(puck, puck.addr, instrs)
  
  weights.puck <- allweights[rownames(allweights) %in% ssg.puck@cell_ID, ]
  
  # figuring out what beads went through l2.rctd 
  l2Ids <- rownames(weights.puck)[(rowSums(subset(weights.puck, select = grep("^l2.", colnames(weights.puck))), na.rm = T) > 0)]
  
  # figuring out which beads only went through l1.rctd
  l1Ids <- rownames(weights.puck)[!rownames(weights.puck) %in% l2Ids]
  
  # Figuring out what cells did not go through l1.rctd
  initFiltIds <- ssg.puck@cell_ID[!ssg.puck@cell_ID %in% rownames(weights.puck)]
  
  # adding the state of the beads to the metadata
  bead.stat <- data.frame(matrix(F, nrow = length(ssg.puck@cell_ID), ncol = 3, dimnames = list(ssg.puck@cell_ID, c('initially.filtered', 'has.l1.rctd', 'has.l2.rctd'))))
  bead.stat[initFiltIds, 'initially.filtered'] <- T
  bead.stat[c(l1Ids, l2Ids), 'has.l1.rctd'] <- T
  bead.stat[l2Ids, 'has.l2.rctd'] <- T
  ssg.puck <- addCellMetadata(ssg.puck, bead.stat)
  
  # adding the weights to the metadata
  w.temp <- data.frame(matrix(nrow = length(ssg.puck@cell_ID), ncol = ncol(weights.puck), 
                              dimnames = list(ssg.puck@cell_ID, colnames(weights.puck))))
  colnames(w.temp) <- colnames(weights.puck)
  w.temp[rownames(weights.puck), colnames(weights.puck)] <- weights.puck
  w.temp[is.na(w.temp)] <- 0
  ssg.puck <- addCellMetadata(ssg.puck, new_metadata = w.temp)
  
  # adding max weight for l1.rctd
  classes <- colnames(ssg.puck@cell_metadata)[grep('^l1.', colnames(ssg.puck@cell_metadata))]
  sel.weights <- subset(ssg.puck@cell_metadata, select = classes)
  maxWeight <- apply(sel.weights, 1, max)
  maxWeight.ct <- colnames(sel.weights)[apply(sel.weights, 1, which.max)]
  maxWeight.ct[ssg.puck@cell_metadata$initially.filtered]
  ssg.puck <- addCellMetadata(ssg.puck, maxWeight, vector_name = 'maxWeight.l1')
  ssg.puck <- addCellMetadata(ssg.puck, maxWeight.ct, vector_name = 'maxCelltype.l1')
  
  # adding max weight for l2.rctd
  classes <- colnames(ssg.puck@cell_metadata)[grep('^l2.', colnames(ssg.puck@cell_metadata))]
  sel.weights <- subset(ssg.puck@cell_metadata, select = classes)
  maxWeight <- apply(sel.weights, 1, max)
  maxWeight.ct <- colnames(sel.weights)[apply(sel.weights, 1, which.max)]
  maxWeight.ct[ssg.puck@cell_metadata$initially.filtered]
  ssg.puck <- addCellMetadata(ssg.puck, maxWeight, vector_name = 'maxWeight.l2')
  ssg.puck <- addCellMetadata(ssg.puck, maxWeight.ct, vector_name = 'maxCelltype.l2')
  
  return(ssg.puck)  
}

mySpatPlot1 <- function(gObj, point_shape = 'no_border', cell_color = NULL, 
                        point_size = 0.4, point_border_stroke = .02, point_alpha = 1, 
                        show_plot=F,save_plot=F, cow_n_col = 1, show_legend=F, cell_color_code = NULL,
                        ctheme = NULL, select_cells = NULL, legend_symbol_size = 4, ...)
{
  plt <- spatPlot(gObj, point_shape=point_shape, cell_color=cell_color, point_size=point_size, 
                  point_border_stroke=point_border_stroke, point_alpha=point_alpha, 
                  show_plot=show_plot, save_plot=save_plot, cow_n_col=cow_n_col, show_legend=show_legend, 
                  cell_color_code=cell_color_code, select_cells=select_cells, 
                  legend_symbol_size=legend_symbol_size, ...) + 
        coord_fixed() 
  
  if (!is.null(ctheme))
    plt <- plt + ctheme
  return(plt)
}

myDimPlot1 <- function(gObj, point_shape = 'no_border', cell_color = NULL, 
                        point_size = 0.4, point_border_stroke = .02, point_alpha = 1, 
                        show_plot=F, show_center_label = F, save_plot=F, cow_n_col = 1, show_legend=T, cell_color_code = NULL,
                        ctheme = NULL, legend_symbol_size = 6, select_cells = NULL, ...)
{
  plt <- dimPlot(gObj, point_shape=point_shape, cell_color=cell_color, point_size=point_size, 
                  point_border_stroke=point_border_stroke, point_alpha=point_alpha, 
                  show_plot=show_plot, show_center_label=show_center_label, save_plot=save_plot, cow_n_col=cow_n_col, show_legend=show_legend, 
                  cell_color_code=cell_color_code, legend_symbol_size=legend_symbol_size, 
                 select_cells=select_cells, legend_symbol_size=legend_symbol_size, ...) + 
    coord_fixed() 
  
  if (!is.null(ctheme))
    plt <- plt + ctheme
  return(plt)
}

theme_noAxis <- function(show_legend = F){
  noAxisTheme <- theme(axis.line=element_blank(),
                       axis.text.x=element_blank(),
                       axis.text.y=element_blank(),
                       axis.ticks=element_blank(),
                       axis.title.x=element_blank(),
                       axis.title.y=element_blank(),
                       # panel.background=element_blank(),
                       panel.border=element_blank(),
                       # panel.grid.major=element_blank(),
                       # panel.grid.minor=element_blank(),
                       # plot.background=element_blank()
  )
  if(!show_legend)
    noAxisTheme <- noAxisTheme + theme(legend.background = element_blank())
  
  return(noAxisTheme)
}

mySpatGenePlot <- function(gObj, gene, slot = 'raw', min.cutoff = NULL, max.cutoff = NULL,
                            point_size = 0.7, axis_title = 8, axis_text = 8, 
                            legend_text = 8, background_color = "white",
                            gradientn_colors = c('gray97', 'orange1', 'orangered4'), 
                            gradientn_vals = c(0, 0.5, 1), hideAxis = T){
  if(slot == 'raw')
    exprs <- gObj@raw_exprs[gene, ]
  else if (slot == 'custom')
    exprs <- gObj@custom_expr[gene, ]
  else if (slot == 'norm')
    exprs <- gObj@norm_expr[gene, ]
  else if (slot == 'scaled')
    exprs <- gObj@norm_scaled_expr[gene, ]
  else
    stop('Wrong slot for gene expression')
  df <- data.frame(x = gObj@spatial_locs$sdimx, y = gObj@spatial_locs$sdimy, expr = exprs)
  df <- df[order(df$expr), ]
  
  if(!is.null(min.cutoff))
    df$expr[df$expr <= min.cutoff] <- min.cutoff
  if(!is.null(max.cutoff))
    df$expr[df$expr >= max.cutoff] <- max.cutoff
  
  plt <- ggplot(df, aes(x, y, color = expr)) + 
    geom_point(size = point_size) + 
    scale_color_gradientn(colors = gradientn_colors, values = gradientn_vals) +
    ggtitle(gene) + 
    coord_fixed() + 
    theme(plot.title = element_text(hjust = 0.5, face = 'bold'), 
          legend.title = element_blank(), legend.text = element_text(size = legend_text), 
          axis.title = element_text(size = axis_title), axis.text = element_text(size = axis_text), 
          panel.grid = element_blank(), panel.background = element_rect(fill = background_color)) 
  if(hideAxis){
    plt <- plt + theme(axis.ticks=element_blank(), axis.text.x=element_blank(),
                       axis.text.y=element_blank(), axis.title = element_blank())
  }
  return(plt)
}


plotGenesRGB3 <- function(gObj, genes, min.expr = NULL, max.expr = NULL, slot = 'raw', 
                          pt.size = 0.5, hideAxis = T, legend = c('bottom', 'right', 'none')[1], 
                          color.hues = c(0.5, 0.83, 0.16), color.values = c(1, 1, 1)){
  # If min.expr is null, then it will be set to the minimum level of expression according to `slot`
  # If min.expr is a vector (length same as genes), expression levels smaller than min.expr are set to min.expr. 
  # Same rule in the opposite direction applies to max.expr
  if(slot == 'raw')
    exprs <- gObj@raw_exprs[genes, ]
  else if (slot == 'custom')
    exprs <- gObj@custom_expr[genes, ]
  else if (slot == 'norm')
    exprs <- gObj@norm_expr[genes, ]
  else if (slot == 'scaled')
    exprs <- gObj@norm_scaled_expr[genes, ]
  else
    stop('Wrong slot for gene expression')
  
  pal1 <- function(x) hsv(color.hues[1], x, color.values[1]) # 1-sqrt(1-x)
  pal2 <- function(x) hsv(color.hues[2], x, color.values[2]) # 1-sqrt(1-x)
  pal3 <- function(x) hsv(color.hues[3], x, color.values[3]) # 1-sqrt(1-x)
  colorMixer <- function(sat.df){
    # columns should be normalized saturation values between 0 and 1
    col1 <- col2rgb(pal1(sat.df[[1]]))
    col2 <- col2rgb(pal2(sat.df[[2]]))
    if (ncol(sat.df) == 2)
      col3 <- col2rgb(pal3(rep_len(0, dim(sat.df)[1]))) # white
    else
      col3 <- col2rgb(pal3(sat.df[[3]]))
    # mixed <- col1 + col2 + col3
    mixed <- 255 - ((255 - col1) + (255 - col2) + (255 - col3))
    mixed[mixed > 255] <- 255
    mixed[mixed < 0] <- 0
    mixed <- mixed / 255
    return(rgb(mixed['red', ], mixed['green', ], mixed['blue', ]))
  }
  plotColorbar <- function(pal, min_val, max_val, orientation = 'horizontal'){
    red <- data.frame(r = seq(min_val, max_val, length.out = 50))
    if(orientation == 'horizontal')
      cbar <- ggplot(red) + geom_tile(aes(x=r, y=0, fill = r)) + 
        scale_fill_gradientn(colours = pal(seq(0, 1, length.out = 50)), guide='none') + 
        theme_bw() + 
        theme(panel.grid = element_blank(), panel.border = element_blank(), 
              axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
              axis.title.y = element_blank()) 
    if(orientation == 'vertical')
      cbar <- ggplot(red) + geom_tile(aes(x=0, y=r, fill = r)) + 
        scale_fill_gradientn(colours = pal(seq(0, 1, length.out = 50)), guide='none') + 
        theme_bw() + 
        theme(panel.grid = element_blank(), panel.border = element_blank(), 
              axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
              axis.title.y = element_blank()) 
    return(cbar)
  }
  
  ng <- length(genes)
  expr.mat <- as.data.frame(as.matrix(t(exprs[genes, ])))
  
  if(is.null(min.expr))
    min.expr <- apply(expr.mat, 2, min)
  if(is.null(max.expr))
    max.expr <- apply(expr.mat, 2, max)
  
  for(i in 1:length(genes)){
    expr.mat[, i][expr.mat[, i] <= min.expr[i]] <- min.expr[i]
    expr.mat[, i][expr.mat[, i] >= max.expr[i]] <- max.expr[i]
    expr.mat[, i] <- (expr.mat[, i] - min.expr[i]) / (max.expr[i] - min.expr[i])
  }
  
  gobj.f <- gObj
  cols <- colorMixer(expr.mat)
  toPlot.df <- cbind(gobj.f@spatial_locs, col=cols)
  toPlot.df <- toPlot.df[order(rowSums(expr.mat)), ]
  
  spPlot <- ggplot() + 
    geom_point(data = gObj@spatial_locs, aes(x = sdimx, y = sdimy), color = 'white', size = .3) + 
    geom_point(data = toPlot.df, aes(x = sdimx, y = sdimy, color = cell_ID), size = pt.size) +
    scale_color_manual(values = setNames(toPlot.df$col, as.character(toPlot.df$cell_ID))) +
    theme_bw() + theme(legend.position = 'none', panel.grid = element_blank()) + coord_fixed()
  
  if(hideAxis){
    spPlot <- spPlot + theme(axis.ticks=element_blank(), axis.text.x=element_blank(),
                             axis.text.y=element_blank(), axis.title = element_blank(), panel.border = element_blank())
  }
  
  cbar.orit <- if (legend == 'bottom') 'horizontal' else 'vertical'
  cyan_plt <- plotColorbar(pal1, min.expr[1], max.expr[1], cbar.orit) + xlab(genes[1])
  mag_plt <- plotColorbar(pal2, min.expr[2], max.expr[2], cbar.orit) + xlab(genes[2])
  if(legend == 'none')
    out = spPlot
  else if(ng > 2){
    yel_plt <- plotColorbar(pal3, min.expr[3], max.expr[3], cbar.orit) + xlab(genes[3])
    if (legend == 'bottom'){
      cbars <- plot_spacer() + cyan_plt + mag_plt + yel_plt + plot_spacer() +
        plot_layout(ncol = 5, widths = c(0.5, 2, 2, 2, 0.5))
      out <- wrap_plots(spPlot, cbars, nrow = 2, heights = c(15, 1))
    }
    else{ # legend == 'right'
      cbars <- plot_spacer() + cyan_plt + mag_plt + yel_plt + plot_spacer() +
        plot_layout(nrow = 5, heights = c(0.5, 1, 1, 1, 0.5))
      out <- wrap_plots(spPlot, cbars, ncol = 2, widths = c(12, 1))
    }
    
  }else{ # two genes/ celltypes
    if (legend == 'right'){
      cbars <- plot_spacer() + cyan_plt + mag_plt + plot_spacer() +
        plot_layout(nrow = 4, heights = c(0.5, 2, 2, 0.5))
      out <- wrap_plots(spPlot, cbars, ncol = 2, widths = c(12, 1))
    }
    else{ # legend == 'bottom'
      cbars <- plot_spacer() + cyan_plt + mag_plt + plot_spacer() +
        plot_layout(ncol = 4, widths = c(0.5, 1, 1, 0.5))
      out <- wrap_plots(spPlot, cbars, nrow = 2, heights = c(15, 1))
    }
  }
  return(out)
}

plotMetadataRGB <- function(gObj, metaCols, min.expr = NULL, max.expr = NULL, 
                            pt.size = 0.5, hideAxis = T, legend = c('bottom', 'right', 'none')[1], 
                            color.hues = c(0.5, 0.83, 0.16), color.values = c(1, 1, 1)){
  # If min.expr is null, then it will be set to the minimum level of expression according to `slot`
  # If min.expr is a vector (length same as metadata), expression levels smaller than min.expr are set to min.expr. 
  # Same rule in the opposite direction applies to max.expr
  exprs <- subset(gObj@cell_metadata, select = metaCols)
  pal1 <- function(x) hsv(color.hues[1], x, color.values[1]) # 1-sqrt(1-x)
  pal2 <- function(x) hsv(color.hues[2], x, color.values[2]) # 1-sqrt(1-x)
  if(length(metaCols) > 2)
    pal3 <- function(x) hsv(color.hues[3], x, color.values[3]) # 1-sqrt(1-x)
  
  colorMixer <- function(sat.df){
    # columns should be normalized saturation values between 0 and 1
    col1 <- col2rgb(pal1(sat.df[[1]]))
    col2 <- col2rgb(pal2(sat.df[[2]]))
    if (ncol(sat.df) == 2)
      col3 <- col2rgb(rep_len("#FFFFFF", dim(sat.df)[1])) # white
    else
      col3 <- col2rgb(pal3(sat.df[[3]]))
    # mixed <- col1 + col2 + col3
    mixed <- 255 - ((255 - col1) + (255 - col2) + (255 - col3))
    mixed[mixed > 255] <- 255
    mixed[mixed < 0] <- 0
    mixed <- mixed / 255
    return(rgb(mixed['red', ], mixed['green', ], mixed['blue', ]))
  }
  plotColorbar <- function(pal, min_val, max_val, orientation = 'horizontal'){
    red <- data.frame(r = seq(min_val, max_val, length.out = 50))
    if(orientation == 'horizontal')
      cbar <- ggplot(red) + geom_tile(aes(x=r, y=0, fill = r)) + 
        scale_fill_gradientn(colours = pal(seq(0, 1, length.out = 50)), guide='none') + 
        theme_bw() + 
        theme(panel.grid = element_blank(), panel.border = element_blank(), 
              axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
              axis.title.y = element_blank()) 
    if(orientation == 'vertical')
      cbar <- ggplot(red) + geom_tile(aes(x=0, y=r, fill = r)) + 
        scale_fill_gradientn(colours = pal(seq(0, 1, length.out = 50)), guide='none') + 
        theme_bw() + 
        theme(panel.grid = element_blank(), panel.border = element_blank(), 
              axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
              axis.title.y = element_blank()) 
    return(cbar)
  }
  
  ng <- length(metaCols)
  expr.mat <- as.data.frame(as.matrix((exprs)))
  
  if(is.null(min.expr))
    min.expr <- apply(expr.mat, 2, min)
  if(is.null(max.expr))
    max.expr <- apply(expr.mat, 2, max)
  
  for(i in 1:length(metaCols)){
    expr.mat[, i][expr.mat[, i] <= min.expr[i]] <- min.expr[i]
    expr.mat[, i][expr.mat[, i] >= max.expr[i]] <- max.expr[i]
    expr.mat[, i] <- (expr.mat[, i] - min.expr[i]) / (max.expr[i] - min.expr[i])
  }
  gobj.f <- gObj
  cols <- colorMixer(expr.mat)
  toPlot.df <- cbind(gobj.f@spatial_locs, col=cols)
  toPlot.df <- toPlot.df[order(rowSums(expr.mat)), ]
  
  spPlot <- ggplot() + 
    geom_point(data = gObj@spatial_locs, aes(x = sdimx, y = sdimy), color = 'white', size = .3) + 
    geom_point(data = toPlot.df, aes(x = sdimx, y = sdimy, color = cell_ID), size = pt.size) +
    scale_color_manual(values = setNames(toPlot.df$col, as.character(toPlot.df$cell_ID))) +
    theme_bw() + theme(legend.position = 'none', panel.grid = element_blank()) + coord_fixed()
  if(hideAxis){
    spPlot <- spPlot + theme(axis.ticks=element_blank(), axis.text.x=element_blank(),
                             axis.text.y=element_blank(), axis.title = element_blank(), panel.border = element_blank())
  }
  cbar.orit <- if (legend == 'bottom') 'horizontal' else 'vertical'
  cyan_plt <- plotColorbar(pal1, min.expr[1], max.expr[1], cbar.orit) + xlab(metaCols[1])
  mag_plt <- plotColorbar(pal2, min.expr[2], max.expr[2], cbar.orit) + xlab(metaCols[2])
  if(legend == 'none')
    out = spPlot
  else if(ng > 2){
    yel_plt <- plotColorbar(pal3, min.expr[3], max.expr[3], cbar.orit) + xlab(metaCols[3])

    if (legend == 'bottom'){
      cbars <- plot_spacer() + cyan_plt + mag_plt + yel_plt + plot_spacer() + 
        plot_layout(ncol = 5, widths = c(0.5, 2, 2, 2, 0.5))
      out <- wrap_plots(spPlot, cbars, nrow = 2, heights = c(15, 1))
    }
    else{ # legend == 'right'
      cbars <- plot_spacer() + cyan_plt + mag_plt + yel_plt + plot_spacer() + 
        plot_layout(nrow = 5, heights = c(0.5, 1, 1, 1, 0.5))
      out <- wrap_plots(spPlot, cbars, ncol = 2, widths = c(12, 1))
    }
    
  }else{ # two genes/ celltypes 
    if (legend == 'right'){
      cbars <- plot_spacer() + cyan_plt + mag_plt + plot_spacer() + 
        plot_layout(nrow = 4, heights = c(0.5, 2, 2, 0.5))
      out <- wrap_plots(spPlot, cbars, ncol = 2, widths = c(12, 1))
    }
    else{ # legend == 'bottom'
      cbars <- plot_spacer() + cyan_plt + mag_plt + plot_spacer() + 
        plot_layout(ncol = 4, widths = c(0.5, 1, 1, 0.5))
      out <- wrap_plots(spPlot, cbars, nrow = 2, heights = c(15, 1))
    }
  }
  return(out)
  
}

getPuckCounts <- function(puck, puck.addr){
  feature.file <- paste(puck.addr, paste0(puck, '.digital_expression_features.tsv'), 
                        sep = '/')
  feature.names <- read.csv(feature.file, header = FALSE, sep = '\t', stringsAsFactors = FALSE)
  
  cell.bc.file <- paste(puck.addr, paste0(puck, '.digital_expression_barcodes.tsv'), 
                        sep = '/')
  cell.bc <- read.csv(cell.bc.file, sep = '\t', header = FALSE, stringsAsFactors = FALSE)$V1
  cell.bc <- sapply(strsplit(cell.bc, split = '-'), function(x) paste0(puck, '-', x[1])) # adding puck name to the bead barcodes
  
  counts <- readMM(paste(puck.addr, paste0(puck, '.digital_expression_matrix.mtx.gz'),
                         sep = '/'))
  dimnames(counts) <- list(feature.names$V2, cell.bc) # rows, columns
  
  return(counts)
}

getPuckLocs <- function(puck, puck.addr){
  bead.coords.file <- paste(puck.addr, 'BeadLocationsForR.csv', sep = '/')
  bead.coords <- read.csv(file = bead.coords.file, header = TRUE, sep = ',')
  return(bead.coords)
}

readSlideSeqGiotto <- function(puck, puck.addr, instrs = NULL, python_path = NULL, plot_dir = './'){
  if(is.null(instrs))
    instrs = createGiottoInstructions(save_dir = plot_dir,
                                      save_plot = FALSE, show_plot = FALSE,
                                      python_path = python_path, plot_format = 'pdf')
  
  counts <- getPuckCounts(puck, puck.addr)
  locs <- getPuckLocs(puck, puck.addr)
  ssg <- createGiottoObject(raw_exprs = counts, spatial_locs = locs, instructions = instrs) # slide-seq Giotto
  return(ssg)
}

subsetRCTDList <- function(rctd.list.results, cellIDs){
  puck.results <- lapply(rctd.list.results, function(res) res[rownames(res) %in% cellIDs, ])
  names(puck.results) <- names(rctd.list.results)
  return(puck.results)
}

addRCTDInfo <- function(gobj, rctd.list, puck, meta.prefix = 'rctd.', addWeights = FALSE){
  puck.cellIDs <- rctd.list$spRNA.cellIDs[grepl(puck, rctd.list$spRNA.cellIDs)]
  puck.results <- subsetRCTDList(rctd.list$results, puck.cellIDs)
  doublets <- puck.results$results_df[, c('spot_class', 'first_type', 'second_type')]
  colnames(doublets) <- paste0(meta.prefix, c('spot_class', 'first_type', 'second_type'))
  ssg.cf <- subsetGiotto(gobj, cell_ids = puck.cellIDs) # filter on counts
  ssg.cf <- addCellMetadata(ssg.cf, doublets)  
  
  if(addWeights){
    weights <- puck.results$weights
    weights <- weights / rowSums(weights) * 100
    colnames(weights) <- paste0(meta.prefix, colnames(weights))
    ssg.cf <- addCellMetadata(ssg.cf, weights)
  }
  return(ssg.cf)
}

addRCTDMaxWeight <- function(gobj, rctd.list, meta.prefix = 'rctd.'){
  rctd.weights <- subsetRCTDList(rctd.list$results, gobj@cell_ID)$weights
  rctd.weights <- rctd.weights / rowSums(rctd.weights)
  maxWeight <- 100 * apply(rctd.weights, 1, max)
  maxWeight.ct <- colnames(rctd.weights)[apply(rctd.weights, 1, which.max)]
  names(maxWeight.ct) <- rownames(rctd.weights)
  
  gobj <- addCellMetadata(gobj, maxWeight.ct, vector_name = paste0(meta.prefix, 'maxCelltype'))
  gobj <- addCellMetadata(gobj, maxWeight, vector_name = paste0(meta.prefix, 'maxWeight'))
  return(gobj)
}

merge.sparse = function(listMatrixes) {
  # adapted from https://xspdf.com/help/52236148.html
  # takes a list of sparse matrixes with different columns and adds them row wise
  
  allColnames <- sort(unique(unlist(lapply(listMatrixes,colnames))))
  for (currentMatrix in listMatrixes) {
    newColLocations <- match(colnames(currentMatrix),allColnames)
    indexes <- which(currentMatrix>0, arr.ind = T)
    newColumns <- newColLocations[indexes[,2]]
    rows <- indexes[,1]
    newMatrix <- sparseMatrix(i=rows,j=newColumns, x=currentMatrix@x,
                              dims=c(max(rows),length(allColnames)))
    rownames(newMatrix) <- rownames(currentMatrix)
    if (!exists("matrixToReturn")) {
      matrixToReturn <- newMatrix
    }
    else {
      matrixToReturn <- rbind2(matrixToReturn,newMatrix)
    }
  }
  colnames(matrixToReturn) <- allColnames
  matrixToReturn  
}
