require(Giotto)
require(ggplot2)
require(Matrix)
require(patchwork)

# A function for adding RCTD weights without removing beads 
readPuck <- function(puck, puck.addr, allweights, instrs, version){
  ssg.puck <- readSlideSeqGiotto2(puck, puck.addr, instrs, version = version)
  
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

# A function for adding RCTD weights without removing beads - for KPMP data
readPuck_j <- function(puck, puck.addr, allweights, instrs){
  ssg.puck <- readSlideSeqGiotto3(puck, puck.addr, instrs)
  
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

filterWeights <- function(gObj, w_column="maxWeight.l2", ct_column="maxCelltype.l2",
                          w_cutoff=50, ct_exception="l2.REN", w_exception=40){
  select_inds <- (gObj@cell_metadata[[w_column]] >= w_cutoff) | 
    ((gObj@cell_metadata[[w_column]] >= w_exception) & 
       gObj@cell_metadata[[ct_column]] %in% ct_exception)
  return(subsetGiotto(gObj, cell_ids = gObj@cell_ID[select_inds]))
}

IsBetweenLines <- function(gObj, lineList){
  # given two lines in lineList and giotto object, returns a boolean vector whether
  # the beads fall between the two lines. Line names have to be "topLine" and "botLine".
  # Each line is specified by two points, p1 and p2
  isBet <- rep_len(FALSE, length(gObj@cell_ID))
  locs.df = gObj@spatial_locs
  if('topLine' %in% names(lineList)){
    topLine = lineList[['topLine']]
    m = (topLine$p2['y'] - topLine$p1['y']) / (topLine$p2['x'] - topLine$p1['x'])
    # y = m * (x - x1) + y1
    locs.df[['topLine']] <- locs.df$sdimy - (m * (locs.df$sdimx - topLine$p1['x']) + topLine$p1['y'])
    locs.df[['isUnderTopLine']] <- locs.df[['topLine']] <= 0
  } else
    locs.df[['isUnderTopLine']] <- TRUE
  if('botLine' %in% names(lineList)){
    botLine = lineList[['botLine']]
    m = (botLine$p2['y'] - botLine$p1['y']) / (botLine$p2['x'] - botLine$p1['x'])
    # y = m * (x - x1) + y1
    locs.df[['botLine']] <- locs.df$sdimy - (m * (locs.df$sdimx - botLine$p1['x']) + botLine$p1['y'])
    locs.df[['isAboveBotLine']] <- locs.df[['botLine']] >= 0
  }
  else
    locs.df[['isAboveBotLine']] <- TRUE
  isBet <- locs.df[['isAboveBotLine']] & locs.df[['isUnderTopLine']]
  return(isBet)
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

plotGenesRGB_DARK <- function(gObj, genes, min.expr = NULL, max.expr = NULL, slot = 'raw', 
                           pt.size = 0.5, hideAxis = T, legend = c('bottom', 'right', 'none')[1]){
  # Dark background. Dots only green, blue or red. 
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
  
  pal1 <- function(x) rgb(x, 0, 0) # hsv(0, 1, 0) # 1-sqrt(1-x)
  pal2 <- function(x) rgb(0, x, 0) # hsv(1/3, 1, 0) # 1-sqrt(1-x)
  pal3 <- function(x) rgb(0, 0, x) # hsv(2/3, 1, 0) # 1-sqrt(1-x)
  
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
  if(length(genes) == 2)
    expr.mat[, 3] <- 0
  
  gobj.f <- gObj
  cols <- rgb(expr.mat)
  toPlot.df <- cbind(gobj.f@spatial_locs, col=cols)
  toPlot.df <- toPlot.df[order(rowSums(expr.mat)), ]
  
  spPlot <- ggplot() + 
    geom_point(data = gObj@spatial_locs, aes(x = sdimx, y = sdimy), color = 'black', size = .3) + 
    geom_point(data = toPlot.df, aes(x = sdimx, y = sdimy, color = cell_ID), size = pt.size) +
    scale_color_manual(values = setNames(toPlot.df$col, as.character(toPlot.df$cell_ID))) +
    theme(legend.position = 'none', panel.grid = element_blank()) + 
    theme(panel.background = element_rect(fill = 'black', color="black")) + 
    coord_fixed()
  
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


plotMetadata_DARK <- function(gObj, metaCols, min.expr = NULL, max.expr = NULL, 
                            pt.size = 0.5, hideAxis = T, legend = c('bottom', 'right', 'none')[1]){
  # If min.expr is null, then it will be set to the minimum level of expression according to `slot`
  # If min.expr is a vector (length same as metadata), expression levels smaller than min.expr are set to min.expr. 
  # Same rule in the opposite direction applies to max.expr
  exprs <- subset(gObj@cell_metadata, select = metaCols)
  pal1 <- function(x) rgb(x, 0, 0) # hsv(0, 1, 0) # 1-sqrt(1-x)
  pal2 <- function(x) rgb(0, x, 0) # hsv(1/3, 1, 0) # 1-sqrt(1-x)
  pal3 <- function(x) rgb(0, 0, x) # hsv(2/3, 1, 0) # 1-sqrt(1-x)
  
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
  if(length(metaCols) == 2)
    expr.mat[, 3] <- 0
  
  gobj.f <- gObj
  cols <- rgb(expr.mat)
  toPlot.df <- cbind(gobj.f@spatial_locs, col=cols)
  toPlot.df <- toPlot.df[order(rowSums(expr.mat)), ]
  
  spPlot <- ggplot() + 
    geom_point(data = gObj@spatial_locs, aes(x = sdimx, y = sdimy), color = 'black', size = .3) + 
    geom_point(data = toPlot.df, aes(x = sdimx, y = sdimy, color = cell_ID), size = pt.size) +
    scale_color_manual(values = setNames(toPlot.df$col, as.character(toPlot.df$cell_ID))) +
    theme(legend.position = 'none', panel.grid = element_blank()) + 
    theme(panel.background = element_rect(fill = 'black', color="black")) + 
    coord_fixed()
  
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

getPuckCounts2 <- function(puck, puck.addr){
  feature.file <- paste(puck.addr, paste0(puck, '.matched.digital_expression_features.tsv.gz'), 
                        sep = '/')
  feature.names <- read.table(feature.file)
  
  cell.bc.file <- paste(puck.addr, paste0(puck, '.matched.digital_expression_barcodes.tsv.gz'), 
                        sep = '/')
  cell.bc <- read.table(cell.bc.file)$V1
  cell.bc <- sapply(strsplit(cell.bc, split = '-'), function(x) paste0(puck, '-', x[1])) # adding puck name to the bead barcodes
  
  counts <- readMM(paste(puck.addr, paste0(puck, '.matched.digital_expression_matrix.mtx.gz'),
                         sep = '/'))
  dimnames(counts) <- list(feature.names$V2, cell.bc) # rows, columns
  
  return(counts)
}

getPuckLocs <- function(puck, puck.addr){
  bead.coords.file <- paste(puck.addr, 'BeadLocationsForR.csv', sep = '/')
  bead.coords <- read.csv(file = bead.coords.file, header = TRUE, sep = ',')
  bead.coords$barcodes <- sapply(bead.coords$barcodes, function(x) paste0(puck, '-', x[1])) # adding puck name to the bead barcodes
  return(bead.coords)
}

getPuckLocs2 <- function(puck, puck.addr){
  bead.coords.file <- paste(puck.addr, paste0(puck, "_barcode_matching.txt.gz"), sep = '/')
  bead.coords <- read.table(file = bead.coords.file, header = F, 
                            col.names = c("sequenced_barcode", "matched_barcode", "xcoord", "ycoord"))
  bead.coords <- unique(bead.coords[, 2:4])  # subset to the matched barcodes
  bead.coords$matched_barcode <- sapply(strsplit(bead.coords$matched_barcode, split = '-'), 
                                        function(x) paste0(puck, '-', x[1])) # adding puck name to the bead barcodes
  return(bead.coords)
}

readSlideSeqGiotto2 <- function(puck, puck.addr, instrs = NULL, python_path = NULL, plot_dir = './', version=1){
  # `version` is the pipeline version. It could be 1 or 2
  if(is.null(instrs))
    instrs = createGiottoInstructions(save_dir = plot_dir,
                                      save_plot = FALSE, show_plot = FALSE,
                                      python_path = python_path, plot_format = 'pdf')
  if(version == 1){
    counts <- getPuckCounts(puck, puck.addr)
    locs <- getPuckLocs(puck, puck.addr)
    rownames(locs) <- locs$barcodes
    locs <- locs[colnames(counts), ] # re-order the barcode locations
    colnames(locs) <- c("barcode", "sdimx", "sdimy")
    locs <- locs[, c("sdimx", "sdimy")]
  }else{
    counts <- getPuckCounts2(puck, puck.addr)
    locs <- getPuckLocs2(puck, puck.addr)
    rownames(locs) <- locs$matched_barcode
    locs <- locs[colnames(counts), ] # re-order the barcode locations
    colnames(locs) <- c("barcode", "sdimx", "sdimy")
    locs <- locs[, c("sdimx", "sdimy")]
  }
  ssg <- createGiottoObject(raw_exprs = counts, spatial_locs = locs, instructions = instrs) # slide-seq Giotto
  return(ssg)
}

readSlideSeqGiotto3 <- function(puck, puck.addr, instrs = NULL, python_path = NULL){
  # function for making Giotto objects for KPMP pucks
  if(is.null(instrs))
    instrs <- createGiottoInstructions(python_path = python_path,
                                      show_plot = F, save_plot = F)

  tab <- read.table(file = paste0(puck.addr, puck, '/MappedDGEForR.csv'), sep = ',', header = T, row.names = 1)
  locs <- read.table(paste0(puck.addr, puck, '/BeadLocationsForR.csv'), sep = ",", header = T)
  locs <- locs[locs$barcodes %in% colnames(tab), ]
  
  # adding puck name to the bead barcode
  colnames(tab) <- paste0(puck, "-", colnames(tab))
  locs$barcodes <- paste0(puck, "-", locs$barcodes)
  
  rownames(locs) <- locs$barcodes
  colnames(locs) <- c("barcode", "sdimx", "sdimy")
  locs <- locs[, c("sdimx", "sdimy")]
  ssg <- createGiottoObject(tab, locs, instructions = instrs)
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
  rctd.weights <- sweep(rctd.weights, 1, rowSums(rctd.weights), '/')
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

stair <- function(x, breaks, levels){
  trans <- rep(0, length(x))
  for(i in seq_along(breaks))
    if(i == 1)
      trans[x < breaks[1]] <- levels[1]
    else
      trans[(x < breaks[i]) & (x >= breaks[i-1])] <- levels[i]
    trans[x >= tail(breaks, 1)] <- tail(levels, 1)
    return(trans)
}

plot_proximity <- function(prox.res, ctCounts, ctCols, min.w = 0.4, max.w = 4, enrich.cutoff = 0.7,
                           break.nct = c(5, 20, 100, 500, 2000), break.size = c(2, 4, 6, 8, 10, 12),
                           background = 'white', label.size = 4, edge.alpha = 0.5, allowSingles = T,
                           minCtCount = 0){
  # ctCols is a named vector for cell type colors
  # ctCounts is a names vector containing the cell type counts, reflected in node size
  # min.w and max.w specify the range of thickness of the edges. The min and max enrichments are linearly mapped to this range.
  # enrich.cutoff is the threshold value on the mean enrichment above which the edges are kept
  # break.nct is a vector specifying the intervals for the cell type counts. 
  # break.size specifies the dot size for each interval specified by break.nct
  # background could be 'white' or 'gray'
  # cell types with counts smaller than minCtCount are omitted 
  
  # prepping enrichment scores for plotting
  CPscores <- prox.res
  CPscores[['cell_1']] <- sapply(strsplit(as.character(CPscores$unified_int), '--'), function(x) x[1])
  CPscores[['cell_2']] <- sapply(strsplit(as.character(CPscores$unified_int), '--'), function(x) x[2])
  
  CPscores <- CPscores[CPscores$enrichm >= enrich.cutoff, ]
  CPscores <- CPscores[CPscores$type_int == 'hetero', ]
  
  CPscores[['enrich_trans']] <- min.w + (CPscores$enrichm - min(CPscores$enrichm)) / (max(CPscores$enrichm) - min(CPscores$enrichm)) * (max.w - min.w)
  
  # prepping the vertex attribute dataframe
  nodes <- data.frame(name = names(ctCols), 
                      ctcol = ctCols[names(ctCols)])
  nodes[['n_occur']] <- ctCounts[rownames(nodes)]
  nodes <- nodes[(nodes$n_occur > 0) & !is.na(nodes$n_occur) , ]

  if(!allowSingles)
    nodes <- nodes[nodes$name %in% unique(CPscores$cell_1) | nodes$name %in% unique(CPscores$cell_2), ]
  
  # removing low count cell types 
  ct2remove <- names(ctCounts)[ctCounts < minCtCount]
  nodes <- nodes[!nodes$name %in% ct2remove, ]
  CPscores <- CPscores[(!CPscores$cell_1 %in% ct2remove) & (!CPscores$cell_2 %in% ct2remove), ]
  
  
  # prep the labels for dot sizes
  labels <- paste0('[', minCtCount, ', ', break.nct[1], ')')
  for(i in seq_along(break.nct)){
    if (i == length(break.nct))
      labels <- append(labels, paste0('[', break.nct[i], ', )'))
    else
      labels <- append(labels, paste0('[', break.nct[i], ', ', break.nct[i+1], ')'))
  }
  nodes[['size_trans']] <- stair(nodes[['n_occur']], break.nct, break.size)
  
  # print(head(CPscores, 50))
  # print(nodes)
  # print(unique(c(CPscores$cell_1, CPscores$cell_2)))
  # print(unique(c(CPscores$cell_1, CPscores$cell_2)) %in% )
  igd <- tryCatch(igraph::graph_from_data_frame(d = CPscores[, c("cell_1", "cell_2", "enrichm", "enrich_trans")], 
                                                directed = F, vertices = nodes), 
                  error = function(err) {
                    print(err)
                    print(nodes)
                    allct <- unique(c(CPscores$cell_1, CPscores$cell_2)) 
                    print('All cell types:')
                    print(allct)
                    print('Cell types not in nodes:')
                    print(allct[!allct %in% nodes$name])
                    })
  igd <- igraph::graph_from_data_frame(d = CPscores[, c("cell_1", 
                                                        "cell_2", "enrichm", "enrich_trans")], directed = F, vertices = nodes)
  
  if(background == 'white'){
    plt <- ggnet2(igd, color = "ctcol", label=T, label.size = label.size, 
                  edge.size = 'enrich_trans', edge.alpha = edge.alpha, size = 'size_trans') + 
      scale_size_manual("#beads", values = break.size, breaks = break.size, labels = labels)
  }
  else if(background == 'gray'){
    plt <- ggnet2(igd, color = "ctcol", label=T, label.size = label.size, label.color = 'white',
                  edge.size = 'enrich_trans', edge.alpha = edge.alpha, edge.color = 'gray90',
                  size = 'size_trans') + 
      scale_size_manual("#beads", values = break.size, breaks = break.size, labels = labels) + 
      theme(panel.background = element_rect(fill = "grey25"))
    
  }
  return(plt)
}
