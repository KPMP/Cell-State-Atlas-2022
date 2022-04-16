library(doParallel)
library(foreach)
library(GGally)
library(igraph)
library(dplyr)

registerDoParallel(8)

source("./utils.R")

plot_mean_proximity <- function(proxList, ctCols, ctCounts, min.w = 0.4, max.w = 4, enrich.cutoff = 0.7,
                                break.nct = c(5, 20, 100, 500, 2000), break.size = c(2, 4, 6, 8, 10, 12),
                                background = 'white', label.size = 4, edge.alpha = 0.5, allowSingles = T,
                                minCtCount = 0){
  # min.w and max.w specify the range of thickness of the edges. The min and max enrichments are linearly mapped to this range.
  # enrich.cutoff is the threshold value on the mean enrichment above which the edges are kept
  # break.nct is a vector specifying the intervals for the cell type counts. 
  # break.size specifies the dot size for each interval specified by break.nct
  # background could be 'white' or 'gray'
  # ctCols is a named vector for cell type colors
  
  # prepping the mean enrichment table
  enrichs <- lapply(proxList, function(x) subset(x, select = c('unified_int', 'enrichm')))
  enrich.df <- data.frame()
  for(i in seq_along(proxList)){
    encr <- enrichs[[i]]
    enrch.row <- data.frame(matrix(encr$enrichm, nrow = 1))
    colnames(enrch.row) <- encr$unified_int
    rownames(enrch.row) <- names(proxList)[i]
    enrich.df <- dplyr::bind_rows(enrich.df, enrch.row)
  }
  enrich.df[is.na(enrich.df)] <- 0
  enrich.mean <- data.table::data.table(unified_int=colnames(enrich.df), enrichm = colMeans(enrich.df))
  enrich.mean[['type_int']] <- ifelse(lapply(strsplit(enrich.mean$unified_int, '--'), function(x) x[1] == x[2]), 'homo', 'hetero')

  plot_proximity(enrich.mean, ctCounts=ctCounts, ctCols=ctCols, min.w=min.w, max.w=max.w,
                 enrich.cutoff = enrich.cutoff, break.nct = break.nct, break.size = break.size,
                 background = background, label.size = label.size, edge.alpha = edge.alpha,
                 allowSingles = allowSingles, minCtCount = minCtCount)
}

consolidateCounts <- function(region_counts, region='Cortex'){
  if(region == 'Cortex'){
    region_counts['PC'] <- sum(region_counts[c('C-PC', 'M-PC', 'IMCD')])
    region_counts['TAL'] <- sum(region_counts[c('C-TAL', 'M-TAL')])
    region_counts <- region_counts[!names(region_counts) %in% c('C-PC', 'M-PC', 'IMCD',
                                                     'C-TAL', 'M-TAL')]
  } else if(region == 'Medulla'){
    region_counts['PC'] <- sum(region_counts[c('C-PC', 'M-PC')])
    region_counts['DTL'] <- sum(region_counts[c('DTL1', 'DTL2', 'DTL3')])
    region_counts['TAL'] <- sum(region_counts[c('C-TAL', 'M-TAL')])
    region_counts <- region_counts[!names(region_counts) %in% 
                                     c('C-PC', 'M-PC', 'C-TAL', 'M-TAL', 'DTL1', 'DTL2', 'DTL3')]
  }
  return(region_counts)
}


colorTable <- data.table::fread("./Cluster_Color_Table.tsv", select = 1:20)
python_path <- "/media/Home_Raid1_Voyager/kian/anaconda3/envs/giotto_200703/bin/python"

gioObjFolder <- c("./deconv_results/GiottoObjects/")

# Plotting
plot.dir <- paste0("./Plots/")
if (!dir.exists(plot.dir))
  dir.create(plot.dir)

plotf <- paste0(plot.dir, 'proxEnrich_perPuck/')
if(!dir.exists(plotf))
  dir.create(plotf)

l1.l2 <- read.table('./deconv_results/HuBMAP_l1-l2_subclasses.tsv', header = T)

# prepping puck names and addressed
pucks <- list.files(gioObjFolder)
pucks <- pucks[startsWith(pucks, "Puck_") & endsWith(pucks, ".rds")]
pucks <- substr(pucks, 1, nchar(pucks) - 4)

# setting up the color table
l1.colors.df <- colorTable[, c('subclass.l1_label', 'subclass.l1_color')]
l1.colors.df <- l1.colors.df[!duplicated(l1.colors.df), ]
l1.colors.df <- l1.colors.df[l1.colors.df$subclass.l1_label %in% l1.l2$name.l1]
l1.colors <- l1.colors.df$subclass.l1_color
names(l1.colors) <- paste0('l1.', l1.colors.df$subclass.l1_label)

l2.colors.df <- colorTable[, c('subclass.l2_label', 'subclass.l2_color')]
l2.colors.df <- l2.colors.df[!duplicated(l2.colors.df), ]
l2.colors.df <- l2.colors.df[l2.colors.df$subclass.l2_label %in% l1.l2$name.l2]
l2.colors <- l2.colors.df$subclass.l2_color
names(l2.colors) <- l2.colors.df$subclass.l2_label

# making color vectors for cortex and medulla based on the cell type merging that we do
colors.ctx <- l2.colors 
colors.ctx <- colors.ctx[!names(colors.ctx) %in% c('C-PC', 'M-PC', 'IMCD',
                                                   'C-TAL', 'M-TAL')]
colors.ctx['PC'] <- l1.colors['l1.PC']
colors.ctx['TAL'] <- l1.colors['l1.TAL']
colors.mdl <- l2.colors 
colors.mdl <- colors.mdl[!names(colors.mdl) %in% c('C-PC', 'M-PC', 
                                                   'C-TAL', 'M-TAL')]
colors.mdl['PC'] <- l1.colors['l1.PC']
colors.mdl['TAL'] <- l1.colors['l1.TAL']

# write.table(rbind(data.frame(ct = names(colors.mdl), color=colors.mdl, region='medulla'),
#                   data.frame(ct = names(colors.ctx), color=colors.ctx, region='cortex')), 
#             file="./proxEnrich_colors.tsv", sep="\t", row.names = F, quote=F)

# Reading the region information
region_info <- readxl::read_xlsx("./SlideSeq_puck_regions.xlsx", )

# Giotto's stuff
ginstuc <- createGiottoInstructions(python_path = python_path,
                                    show_plot = F, save_plot = F)


prox.list <- foreach(puck=unlist(pucks), .final=function(x) setNames(x, unlist(pucks))) %dopar%{
  reg <- region_info$`Region Level 2`[region_info$Specimen == puck]
 
  ssg <- readRDS(paste0(gioObjFolder, puck, ".rds"))
  
  # limiting to beads with confident l2 calls 
  ssg <- filterWeights(ssg, 'maxWeight.l2', 'maxCelltype.l2', 30, c())
  ssg <- subsetGiotto(ssg, cell_ids = ssg@cell_ID[ssg@cell_metadata$inBoundary])
  
  # removing beads that change cell type between l1 and l2 deconvs
  ssg <- subsetGiotto(ssg, cell_ids = ssg@cell_ID[substr(ssg@cell_metadata$maxCelltype.l2, 4, 100) %in% l1.l2$name.l2])
  
  # defining a modified metadata column for proximity enrichments and a corresponding color vector
  prox.meta <- substr(ssg@cell_metadata$maxCelltype.l2, 4, 100)
  if(reg == "Cortex"){
    prox.meta[prox.meta %in% c('C-PC', 'M-PC', 'IMCD')] <- 'PC'
    prox.meta[prox.meta %in% c('C-TAL', 'M-TAL')] <- 'TAL'
    prox.meta[prox.meta == 'EC-DVR'] <- 'EC-AEA'
    prox.meta[prox.meta  %in% c('EC-AVR', 'ATL', 'DTL1', 'DTL2', 'DTL3', 'M-IC-A', "M-FIB")] <- NA
    prox.meta[prox.meta %in% c('cycPT', "cycCNT", "cycDCT", 'cycEC', 'cycMYOF')] <- NA
    color.vec <- colors.ctx
    }
  else if (reg %in% c("Outer Medulla", "Inner/Outer Medulla", "Inner Medulla")){
    prox.meta[prox.meta %in% c('C-PC', 'M-PC')] <- 'PC'
    prox.meta[prox.meta %in% c('C-TAL', 'M-TAL')] <- 'TAL'
    prox.meta[prox.meta %in% c('DTL1', 'DTL2', "DTL3")] <- 'DTL'
    prox.meta[prox.meta == 'EC-AEA'] <- 'EC-DVR'
    prox.meta[prox.meta  %in% c('DCT', 'PT-S1', 'PT-S2', 'PT-S3', "CNT",
                                'aPT', "POD", "EC-PTC", "EC-GC", "C-IC-A", 'PEC',
                                "MC")] <- NA
    prox.meta[prox.meta %in% c('cycPT', "cycCNT", "cycDCT", 'cycEC', 'cycMYOF')] <- NA
    color.vec <- colors.mdl
  }else
    color.vec <- l2.colors
  
  ssg <- addCellMetadata(ssg, prox.meta, "prox.meta")
  ssg <- subsetGiotto(ssg, ssg@cell_ID[!is.na(ssg@cell_metadata$prox.meta)])
  
  ssg <- createSpatialNetwork(ssg, method = 'Delaunay', minimum_k = 1, maximum_distance_delaunay = 100)
  print(paste0(puck, " started"))
  prox.enrc <- cellProximityEnrichment(ssg, cluster_column = 'prox.meta',
                                       adjust_method = 'fdr', number_of_simulations = 2500)$enrichm_res
  
  n_interacts <- nrow(prox.enrc[original >= 5 & simulations >= 0.5 & 
                                  (p_higher_orig <= 0.05 | p_lower_orig <= 0.05)])
  # barplot
  h = as.integer(16 / 116 * n_interacts) # experimental height value
  pdf(file = paste0(plotf, 'prox_barplot-', puck, '.pdf'), height = h, width = 12)
  print(cellProximityBarplot(gobject = ssg, list(enrichm_res=prox.enrc),
                             show_plot = F, save_plot = F, min_sim_ints = .5))
  dev.off()
  
  # proximity network
  ct.counts <- table(ssg@cell_metadata[['prox.meta']])
  pdf(file = paste0(plotf, 'prox_network-', puck, '.pdf'))
  print(plot_proximity(prox.enrc, ctCounts = ct.counts, 
                       ctCols = color.vec, min.w = 0.4,
                       enrich.cutoff = 1.0, allowSingles = F))
  dev.off()
  
  rm(ssg)
  gc()
  prox.enrc
}
saveRDS(prox.list, "proxEnrich_results.rds")


# calculating cell type counts for pucks
ct.counts <- foreach(puck=pucks, .final=function(x) setNames(x, pucks)) %dopar%{
  reg <- region_info$`Region Level 2`[region_info$Specimen == puck]
  ssg <- readRDS(paste0(gioObjFolder, puck, ".rds"))
  ssg <- filterWeights(ssg, 'maxWeight.l2', 'maxCelltype.l2', 30, c())
  ssg <- subsetGiotto(ssg, cell_ids = ssg@cell_ID[ssg@cell_metadata$inBoundary])
  ssg <- subsetGiotto(ssg, cell_ids = ssg@cell_ID[substr(ssg@cell_metadata$maxCelltype.l2, 4, 100) %in% l1.l2$name.l2])
  df <- data.frame(table(ssg@cell_metadata$maxCelltype.l2))
  colnames(df) <- c('celltype', puck)
  df
}
ctCounts.df <- Reduce(function(x, y) full_join(x, y, by = 'celltype'), ct.counts)
ctCounts.df[is.na(ctCounts.df)] <- 0



reg <- "Cortex"

reg_pucks <- region_info$Specimen[region_info$`Region Level 1`==reg]
reg_counts <- setNames(rowSums(ctCounts.df[, reg_pucks]), substr(ctCounts.df$celltype, 4, 100))
reg_counts <- consolidateCounts(reg_counts, reg)

pdf(paste0(plot.dir, "ED5_d.pdf"), height = 7/10 * 9, width = 9)
plot_mean_proximity(prox.list[reg_pucks], ctCols = colors.ctx, ctCounts = reg_counts, allowSingles = F, 
                    enrich.cutoff = 0.7, min.w=0.25, minCtCount = 100,
                    break.nct = c(200, 500, 2000, 5000, 15000, 25000), 
                    break.size = c(2, 4, 6, 8, 10, 12, 14))
dev.off()


reg <- "Medulla"

reg_pucks <- region_info$Specimen[region_info$`Region Level 1`==reg]
reg_counts <- setNames(rowSums(ctCounts.df[, reg_pucks]), substr(ctCounts.df$celltype, 4, 100))
reg_counts <- consolidateCounts(reg_counts, reg)

pdf(paste0(plot.dir, "ED5_e.pdf"), height = 7/10 * 9, width = 9)
plot_mean_proximity(prox.list[reg_pucks], ctCols = colors.mdl, ctCounts = reg_counts, allowSingles = F, 
                    enrich.cutoff = 0.8, min.w=0.25, minCtCount = 200,
                    break.nct = c(500, 2000, 5000, 10000, 25000, 50000), 
                    break.size = c(2, 4, 6, 8, 10, 12, 14))
dev.off()