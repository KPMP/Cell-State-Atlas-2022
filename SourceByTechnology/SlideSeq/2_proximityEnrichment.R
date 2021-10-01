library(ggplot2)
library(Matrix)
library(Giotto)
library(reshape2)
library(patchwork)
source('utils.R')

nSignificantInteractions <- function(CPscore, min_orig_ints = 5, min_sim_ints = 5, p_val = 0.05){
  table_mean_results_dc = CPscore$enrichm_res
  original = simulations = p_higher_orig = p_lower_orig = enrichm = type_int = unified_int = NULL
  table_mean_results_dc_filter = table_mean_results_dc[original >= 
                                                         min_orig_ints & simulations >= min_sim_ints]
  table_mean_results_dc_filter = table_mean_results_dc_filter[p_higher_orig <= 
                                                                p_val | p_lower_orig <= p_val]
  return(nrow(table_mean_results_dc_filter))
}

python_path <- "/media/Home_Raid1_Voyager/kian/anaconda3/envs/giotto_200703/bin/python"
out.dir <- "./outputs"
slideseq.dir <- "../ForHubmapCollaborators"
rctd.dir <- paste0(out.dir, "/", "deconv_files")
plot.dir <- paste0(out.dir, "/", "Plots")

pucks <- list.dirs(slideseq.dir, full.names = FALSE)
pucks <- pucks[nchar(pucks) > 0]

sel.pucks <- c("Puck_200903_26", "Puck_200903_02", "Puck_200903_01", "Puck_200903_03",
               "Puck_200903_07", "Puck_200903_05", "Puck_200906_18", "Puck_200906_10",
               "Puck_200903_06", "Puck_200903_25")


weights.all <- read.table(paste0(rctd.dir, "/", "combinedWeights.tsv"),
                          sep = '\t', header = T, check.names = F)

instrs = createGiottoInstructions(save_dir = plot.dir,
                                  save_plot = TRUE,
                                  show_plot = FALSE,
                                  python_path = python_path, plot_format = 'pdf')

l1.l2 <- read.table(paste0(rctd.dir, "/", "rctd_subclasses.tsv"), sep = '\t', header = T)

# consolidating the immune cells in subclass l2
imms <- l1.l2$name.l2[l1.l2$name.l1 == 'IMM']
weights.all[, paste0('l2.', imms)] <- NULL
weights.all$l2.IMM <- weights.all$l1.IMM

prox.list <- list()
for(puck in sel.pucks){
  print(puck)
  ssg.puck <- readPuck(puck, paste(slideseq.dir, puck, sep = '/'), weights.all, instrs)
  ssg.puck@cell_metadata$maxCelltype.l2 <- substr(ssg.puck@cell_metadata$maxCelltype.l2, 4, 100)
  ssg.puck.f <- subsetGiotto(ssg.puck, cell_ids = ssg.puck@cell_ID[ssg.puck@cell_metadata$maxWeight.l2 >= 50])
  ssg.puck.f <- createSpatialNetwork(ssg.puck.f, method = 'Delaunay', minimum_k = 1)
  prox.list[[puck]] <- cellProximityEnrichment(ssg.puck.f, cluster_column = 'maxCelltype.l2',
                                               adjust_method = 'fdr', number_of_simulations = 5000)
}

saveRDS(prox.list, file = paste0(out.dir, "/", "goodpucks_prox_l2_imm.rds"))

# Plotting the mean proximity enrichment network ####
prox.list <- readRDS(paste0(out.dir, "/", "goodpucks_prox_l2_imm.rds"))
pucks <- names(prox.list)
sel.pucks <- c("Puck_200903_26", "Puck_200903_02", "Puck_200903_01", "Puck_200903_03", 
               "Puck_200903_07", "Puck_200903_05", "Puck_200906_18", "Puck_200906_10",
               "Puck_200903_06", "Puck_200903_25")

ssg.puck <- readPuck(sel.pucks[1], paste(slideseq.dir, sel.pucks[1], sep = '/'), weights.all, instrs) # dummy Giotto object

enrichs <- lapply(prox.list, function(x) subset(x$enrichm_res,
                                                select = c('unified_int', 'enrichm')))
enrich.df <- data.frame()
for(i in seq_along(pucks)){
  encr <- enrichs[[i]]
  enrch.row <- data.frame(matrix(encr$enrichm, nrow = 1))
  colnames(enrch.row) <- encr$unified_int
  rownames(enrch.row) <- pucks[[i]]
  enrich.df <- dplyr::bind_rows(enrich.df, enrch.row)
}
enrich.df[is.na(enrich.df)] <- 0
enrich.df <- enrich.df[sel.pucks, ]
enrich.mean <- data.table::data.table(unified_int=colnames(enrich.df), enrichm = colMeans(enrich.df))
enrich.mean[['type_int']] <- ifelse(lapply(strsplit(enrich.mean$unified_int, '--'), function(x) x[1] == x[2]), 'homo', 'hetero')
enrich.mean <- enrich.mean[enrich.mean$enrichm > 0.6]
prox.mean <- prox.list[[1]]
prox.mean$enrichm_res <- enrich.mean
pdf(paste0(plot.dir, '/', 'goodpucks_proxNetwork_mean.pdf'))
cellProximityNetwork(ssg.puck, prox.mean, remove_self_edges = T, 
                     layout = 'Fruchterman', edge_width_range = c(1, 3), 
                     node_size = 2, node_text_size = 4,
                     rescale_edge_weights = F, show_plot = F, save_plot = F, 
                     only_show_enrichment_edges = T)
dev.off()