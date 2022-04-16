library(Seurat)
library(stringr)
library(ggplot2)
library(readxl)


all_merged <- readRDS('all_merged_subclass.l2_norazor.RDS')



genes <- c('NPHS2','EMCN','POSTN','VCAM1','SLC5A12','PALMD','MYH11','REN','NOS1')
Idents(all_merged) <- all_merged$transfer_subset
comp <- subset(all_merged,idents = c('POD','EC-GC','MC','PEC','PT-S1','PT-S2',
                                     'EC-AEA','VSMC','MD'))
comp$transfer_subset <- ifelse(comp$transfer_subset == 'PT-S1'|
                                 comp$transfer_subset == 'PT-S2',
                               'PT-S1/2',comp$transfer_subset)
Idents(comp) <- factor(comp$transfer_subset,levels = c('POD','EC-GC','MC','PEC','PT-S1/2',
                                                       'EC-AEA','VSMC','MD'))
DotPlot(comp,'SCT',rev(genes),cols = c('lightgray','purple'),col.min = 0)+
  coord_flip()+
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5))#+
