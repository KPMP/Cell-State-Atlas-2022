library(Seurat)
library(ggplot2)
library(reshape2)
library(readxl)

sample_list <- c('V19S25-019_XY02-M32','V19S25-019_XY03-M61',
                 'V19S25-016_XY01_18-0006','V19S25-019_XY04-F52',
                 'V10S15-102_XY03_IU-21-015-2','V10S15-102_XY02_IU-21-019-5',
                 'V19S25-017_XY03-13437',
                 'V10S14-085_XY01_20-0038','V10S14-085_XY02_20-0040',
                 'V10S14-085_XY03_21-0056','V10S14-085_XY04_21-0057',
                 'V10S14-086_XY01_21-0055','V10S14-086_XY02_20-0039',
                 'V10S14-086_XY03_21-0063','V10S14-086_XY04_21-0066',
                 'V10S14-087_XY01_21-0061','V10S14-087_XY02_21-0063',
                 'V10S14-087_XY03_21-0064','V10S14-087_XY04_21-0065',
                 'V10S21-388_XY01_21-0068','V10S21-388_XY02_20-0071',
                 'V10S21-388_XY03_20-0072','V10S21-388_XY04_20-0073')

for (sample in sample_list){
  spatial <- readRDS(paste0('subclass.l2/',sample,'_seurat_only.RDS'))
  spatial@meta.data$orig.ident <- sample
  DefaultAssay(spatial) <- 'SCT'
  assign(sample, spatial)
}

all_merged <- merge(`V19S25-019_XY02-M32`,`V19S25-019_XY03-M61`)
var_features <- c(VariableFeatures(`V19S25-019_XY02-M32`),VariableFeatures(`V19S25-019_XY03-M61`))
for (sample in sample_list[3:23]){
  all_merged <- merge(all_merged,get(sample))
  var_features <- c(var_features,VariableFeatures(get(sample)))
}

all_merged <- SCTransform(all_merged,assay = 'Spatial')
DefaultAssay(all_merged) <- "SCT"
VariableFeatures(all_merged) <- var_features
all_merged <- RunPCA(all_merged, verbose = FALSE)
all_merged <- FindNeighbors(all_merged, dims = 1:30)
all_merged <- FindClusters(all_merged, verbose = FALSE)
all_merged <- RunUMAP(all_merged, dims = 1:30)

all_merged@meta.data$border <- NA
for (smp in sample_list){
  border <- read.csv(paste0('talmod/',smp,'.border.csv'))
  all_merged@meta.data[all_merged@meta.data$orig.ident == smp,]$border <- as.character(border$border)
}
any(is.na(all_merged@meta.data$border))
Idents(all_merged) <- all_merged@meta.data$border
all_merged <- subset(all_merged,idents = 'razor',invert=T)

akismp <- c('V10S14-085_XY01_20-0038','V10S14-085_XY02_20-0040',
            'V10S14-085_XY03_21-0056',
            'V10S14-086_XY04_21-0066','V10S14-087_XY04_21-0065',
            'V10S21-388_XY02_20-0071')

ckdsmp <- c('V19S25-017_XY03-13437',
            'V10S14-085_XY04_21-0057',
            'V10S14-086_XY01_21-0055','V10S14-086_XY02_20-0039',
            'V10S14-086_XY03_21-0063',
            'V10S14-087_XY01_21-0061','V10S14-087_XY02_21-0063',
            'V10S14-087_XY03_21-0064',
            'V10S21-388_XY01_21-0068',
            'V10S21-388_XY03_20-0072','V10S21-388_XY04_20-0073')

all_merged@meta.data$condition <- ifelse(all_merged@meta.data$orig.ident %in% akismp,'AKI',
                                         ifelse(all_merged@meta.data$orig.ident %in% ckdsmp,
                                                'CKD','Ref'))
all_merged@meta.data$condition.2 <- ifelse(all_merged@meta.data$orig.ident %in% akismp,0,
                                           ifelse(all_merged@meta.data$orig.ident %in% ckdsmp,
                                                  1,0))


saveRDS(all_merged,'all_merged_subclass.l2_norazor.RDS')