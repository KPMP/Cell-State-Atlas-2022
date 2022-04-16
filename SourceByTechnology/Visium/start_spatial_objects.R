#Script to read data and map label from reference dataset

library(Seurat)
library(SeuratDisk)

KBR <- LoadH5Seurat("../WashU-UCSD_SN_Kidney_KPMP-Biopsy_HuBMAP-Ref_10X-R_Seurat_02192021.h5Seurat")

#Update to newest meta data (including cell type annotations)
meta <- read.delim("../WashU-UCSD_SN_Kidney_KPMP-Biopsy_HuBMAP-Ref_10X-R_04162021_Meta.txt",sep="\t",header=TRUE,row.names=1)
KBR@meta.data <- meta[rownames(KBR@meta.data),]
head(KBR@meta.data)

#Pass sample name and label to map via command line
args = commandArgs(trailingOnly=TRUE)
sample <- args[1]
label <- args[2]

Idents(KBR) <- KBR@meta.data[[label]]

spatial <- readRDS(paste0('../../spatial/',sample,'.RDS'))

#Read and normalize spatial data
spatial <- Load10X_Spatial(paste0(sample,'/outs/'))
spatial <- SCTransform(spatial, assay = "Spatial", verbose = FALSE)
spatial <- SCTransform(spatial, assay = "Spatial", verbose = FALSE)
spatial <- RunPCA(spatial, assay = "SCT", verbose = FALSE)
spatial <- FindNeighbors(spatial, dims = 1:30)
spatial <- FindClusters(spatial, verbose = FALSE)
spatial <- RunUMAP(spatial, dims = 1:30)

#Map labels from reference to spots
anchors <- FindTransferAnchors(reference = KBR, query = spatial, normalization.method = "SCT",query.assay='SCT')
predictions.assay <- TransferData(anchorset = anchors, refdata = KBR@meta.data[[label]], prediction.assay = TRUE, 
                                  weight.reduction = spatial[["pca"]])
spatial[["predictions"]] <- predictions.assay

DefaultAssay(spatial) <- "predictions"
df_pred <- spatial@assays[["predictions"]]@data
max_pred <- apply(df_pred,2,function(x) max.col(t(x),'first'))
max_pred_val <- apply(df_pred,2,function(x) max(t(x)))
max_pred <- as.data.frame(max_pred)
max_pred$Seurat_subset <- rownames(df_pred)[max_pred$max_pred]
max_pred$score <- max_pred_val
max_pred$Barcode <- rownames(max_pred)
write.csv(max_pred[,c('Barcode','Seurat_subset')],paste0(sample,'.transfer','.csv'),quote=F,row.names = F)


spatial@meta.data$transfer_subset <- max_pred$Seurat_subset
spatial@meta.data$transfer_subset_score <- max_pred$score

saveRDS(spatial,paste0(sample,'_seurat_only.RDS'))
