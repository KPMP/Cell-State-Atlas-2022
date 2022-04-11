#Script to generate prediction scores feature plots and piecharts plots

library(Matrix)
library(data.table)
library(Seurat)
library(SeuratData)
library(dplyr)
library(gt)
library(igraph)
library(RColorBrewer)
library(cowplot)
library(readxl)
library(stringr)
library(ggplot2)
library(reshape2)

label <- 'subclass.l2'

#Table with celltypes labels and colors

coltable <- read_excel('UCSD_data_032921/Kidney_Reference_Atlas_Summary_Tables.xlsx',
                       sheet='Cluster Color Table',
                       range = 'H1:J101')
coltable <- unique(coltable)

sample <- 'V19S25-017_XY03-13437'
#Size factor varies for each sample
sz=1.5

smp_split <- unlist(str_split(str_remove(sample,'_slide1'),'_'))
img_folder <- smp_split[1]
img_file <- list.files(paste0('../spatial_samples/',img_folder),paste0(smp_split[2],'.*tif'))

spatial <- readRDS(paste0(label,'/',sample,'_seurat_only.RDS'))

dir.create(file.path(label, 'prelim_seurat_only'), showWarnings = FALSE)

getPalette <- colorRampPalette(brewer.pal(9, "YlOrRd"))
cell_types_all <- coltable$subclass.l2_label


spatial@meta.data$transfer_subset <- factor(spatial@meta.data$transfer_subset,
                                            levels = coltable[[1]])
Idents(spatial) <- spatial@meta.data$transfer_subset
pdf(paste0(label,'/prelim_seurat_only/',sample,'_seurat.pdf'),width=12)
plot(SpatialDimPlot(spatial))+
  scale_fill_manual(values=coltable[[3]],limits=coltable[[1]])
dev.off()
preds <- t(as.data.frame(spatial@assays$predictions@data))
colnames(preds) <- paste0('seurat_',colnames(preds))
spatial@meta.data <- cbind(spatial@meta.data,preds)

pdf(paste0(label,'/prelim_seurat_only/',sample,'_features.pdf'),width=8,height = 4)
for (ctype in cell_types_all){
  p1 <- SpatialFeaturePlot(spatial,ctype,crop = T,pt.size.factor = sz)+ 
    scale_fill_gradientn(colors=rev(brewer.pal(n = 11, name = "Spectral")),
                         limits=c(0,0.25),
                         oob = scales::squish)
  plot(p1)
}
dev.off()

#Plot specific gene expression as a feature plot
DefaultAssay(spatial) <- 'SCT'
pdf(paste0(label,'/prelim_seurat_only/',sample,'_SLC12A1.pdf'),width=8,height = 4)
SpatialFeaturePlot(spatial,'SLC12A1',crop = T)+ 
  scale_fill_gradientn(colors=rev(brewer.pal(n = 11, name = "Spectral")),
                       limits=c(1,4.05),
                       oob = scales::squish)
dev.off()


#spript to generate an plot piecharts adapted from SPOTlight source
metadata_ds <- as.data.frame(t(spatial@assays$predictions@data))

rownames(coltable) <- coltable[[1]]
coltable <- coltable[rownames(coltable) %in% colnames(metadata_ds),]
metadata_ds <- metadata_ds[,coltable[[1]]]

rownames(coltable) <- coltable[[1]]
cols <- coltable[colnames(metadata_ds)[colSums(metadata_ds) > 0],][[3]]


spatial_coord <- data.frame(spatial@images[['slice1']]@coordinates) %>%
  tibble::rownames_to_column("barcodeID") %>%
  dplyr::inner_join(metadata_ds %>% tibble::rownames_to_column("barcodeID"),
                    by = "barcodeID")



img <- tiff::readTIFF(paste0('../spatial_samples/',img_folder,'/',img_file))

# Convert image to grob object
img_grob <- grid::rasterGrob(img,
                             interpolate = FALSE,
                             width = grid::unit(1, "npc"),
                             height = grid::unit(1, "npc"))

## Plot spatial scatterpie plot
scatterpie_plt <- suppressMessages(
  ggplot2::ggplot() +
    ggplot2::annotation_custom(
      grob = img_grob,
      xmin = 0,
      xmax = ncol(img),
      ymin = 0,
      ymax = -nrow(img)) +
    scatterpie::geom_scatterpie(
      data = spatial_coord,
      ggplot2::aes(x = imagecol,
                   y = imagerow),
      cols = coltable[[1]],
      color = NA,
      alpha = 1,
      pie_scale = .2*sz) +
    ggplot2::scale_y_reverse() +
    ggplot2::ylim(nrow(img), 0) +
    ggplot2::xlim(0, ncol(img)) +
    cowplot::theme_half_open(11, rel_small = 1) +
    ggplot2::theme_void() +
    ggplot2::coord_fixed(ratio = 1,
                         xlim = NULL,
                         ylim = NULL,
                         expand = TRUE,
                         clip = "on")+
    ggplot2::scale_fill_manual(values = cols))
print(paste0(label,paste0('/prelim_seurat_only/'),sample,'_pie.pdf'))
pdf(paste0(label,paste0('/prelim_seurat_only/'),sample,'_pie.pdf'),width=12,height = 8)
plot(scatterpie_plt)
dev.off()

