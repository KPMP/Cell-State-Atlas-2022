library(Seurat)
library(readxl)
library(reshape2)
library(ggplot2)
library(clustree)
library(stringr)
library(cowplot)

types_table <- read_excel('UCSD_data_032921/Kidney_Reference_Atlas_Summary_Tables.xlsx',
                          sheet='Annotations sn10X',range = 'I3:M103')
types_table <- as.data.frame(unique(types_table))
rownames(types_table) <- types_table$Subclass.l2

  all_merged <- readRDS('all_merged_subclass.l2_norazor.RDS')

preds <- as.data.frame(all_merged@assays$predictions@data)[1:74,]
preds <- t(t(preds)/colSums(preds))
preds <- melt(preds)
preds$subclass.l1 <- types_table[as.character(preds$Var1),'subclass.l1']

preds <- aggregate(preds$value,by=list(Spot=preds$Var2,Subclass.l1=preds$subclass.l1),FUN=sum)
preds <- dcast(preds,Subclass.l1 ~ Spot)
rownames(preds) <- preds$Subclass.l1
preds$Subclass.l1 <- NULL
added_subclass <- rownames(preds)[max.col(t(preds))]
names(added_subclass) <- colnames(preds)

all_merged@meta.data$added_subclass <- added_subclass[rownames(all_merged@meta.data)]
Idents(all_merged) <- all_merged@meta.data$added_subclass



tals <- subset(all_merged,idents = 'TAL')
tals@assays$ctypes <- tals@assays$predictions
DefaultAssay(tals) <- 'ctypes'
tals <- subset(tals,
                     features = row.names(tals@assays$ctypes@data)[1:74])

VariableFeatures(tals) <- row.names(tals@assays$ctypes)
tals <- ScaleData(tals)
tals <- RunPCA(tals)


tals <- FindNeighbors(tals,reduction = 'pca',dims=1:30)
tals <- FindClusters(tals,resolution = 0.45)
tals <- RunUMAP(tals,reduction = 'pca',dims=1:30)

Idents(tals) <- tals@meta.data$ctypes_snn_res.0.45
#saveRDS(tals,'tals_with_umap_after_sw.RDS')
#tals <- readRDS('tals_with_umap_after_sw.RDS')

pdf('niches/umap_tal_sw.pdf')
DimPlot(tals,label = T)
dev.off()
coltable <- read_excel('UCSD_data_032921/Kidney_Reference_Atlas_Summary_Tables.xlsx',
                       sheet='Misc Color Table',
                       range = 'H1:J4')
coltable <- unique(coltable)
Idents(tals) <- factor(tals@meta.data$condition,levels = coltable$condition.l1_label)
pdf('niches/umap_tal_condition.pdf')
DimPlot(tals,label = F,cols = coltable$condition.l1_color)
dev.off()
pdf('niches/umap_tal_patient.pdf',width = 12)
DimPlot(tals,label = F,group.by = 'orig.ident')
dev.off()
pdf('niches/dotplot_tal_sw_all.pdf',width = 12,height = 6)
DotPlot(tals,features = c('C-TAL','M-TAL','aTAL1',
                          'aTAL2','dC-TAL','dM-TAL',
                          types_table[types_table$subclass.l1 == 'FIB',]$Subclass.l2,
                          types_table[types_table$subclass.l1 == 'IMM',]$Subclass.l2),
        col.min = 0)+
  theme(axis.text.x = element_text(angle = 90,vjust = .5))
dev.off()

preds <- as.data.frame(tals@assays$predictions@data)[1:74,]
preds <- t(t(preds)/colSums(preds))
tals[['test']] <- tals[['predictions']]
tals@assays$test@data <- log(preds +1)
pdf('niches/dotplot_zscore_tal_sw_all.pdf',height = 4)
DotPlot(tals,features = c('C-TAL','M-TAL','aTAL1',
                          'aTAL2','dC-TAL','dM-TAL','MD',
                          types_table[types_table$subclass.l1 == 'FIB',]$Subclass.l2,
                          types_table[types_table$subclass.l1 == 'IMM',]$Subclass.l2),
        assay = 'test')+
  theme(axis.text.x = element_text(angle = 90,vjust = .5))
dev.off()


coltable <- read_excel('UCSD_data_032921/Kidney_Reference_Atlas_Summary_Tables.xlsx',
                       sheet='Cluster Color Table',
                       range = 'H1:J101')
coltable <- unique(coltable)
coltable <- coltable[coltable$subclass.l2_label %in%
                       c('aTAL1','aTAL2',
                         'M-TAL','dM-TAL',
                         'C-TAL','dC-TAL','MD',
                         types_table[types_table$subclass.l1 == 'FIB',]$Subclass.l2,
                         types_table[types_table$subclass.l1 == 'IMM',]$Subclass.l2),]


talpreds <- as.data.frame(tals@assays$predictions@data)[1:74,]
talpreds <-t(t(talpreds)/colSums(talpreds))
talpreds <- melt(talpreds)
talpreds <- talpreds[talpreds$Var1 %in%
                     c('aTAL1','aTAL2',
                       'M-TAL','dM-TAL',
                       'C-TAL','dC-TAL','MD',
                       types_table[types_table$subclass.l1 == 'FIB',]$Subclass.l2,
                       types_table[types_table$subclass.l1 == 'IMM',]$Subclass.l2),]
talpreds <- as.data.frame(talpreds)
talpreds$cluster <- tals@meta.data[talpreds$Var2,'ctypes_snn_res.0.45']
talpreds$Proportion <- NA
for (cl in unique(talpreds$cluster)){
  talpreds[talpreds$cluster == cl, 'Proportion'] <- talpreds[talpreds$cluster == cl, 'value'] /
    length(unique(talpreds[talpreds$cluster == cl,'Var2']))
}
talpreds <- aggregate(talpreds$Proportion,by=list(Type=talpreds$Var1,Cluster=talpreds$cluster),FUN=sum)
talpreds <- talpreds[talpreds$x > 0,]


talpreds$Type <- factor(talpreds$Type,levels = coltable$subclass.l2_label)

pdf('niches/barplot_tal_sw_relative_cluster_proportions.pdf')
ggplot(talpreds,aes(x=Cluster,y=x,fill=Type))+
  geom_bar(stat = 'identity',position = 'stack')+
  scale_fill_manual(values = coltable$subclass.l2_color)
dev.off()  


talpreds <- as.data.frame(tals@assays$predictions@data)[1:74,]
talpreds <-t(t(talpreds)/colSums(talpreds))
talpreds <- melt(talpreds)
talgroup <- c('aTAL1','aTAL2',
              'M-TAL','dM-TAL',
              'C-TAL','dC-TAL','MD',
              types_table[types_table$subclass.l1 == 'FIB',]$Subclass.l2,
              types_table[types_table$subclass.l1 == 'IMM',]$Subclass.l2)
talrows <- !talpreds$Var1 %in% talgroup
talpreds$Var1 <- as.character(talpreds$Var1)
talpreds[talrows,'Var1'] <- types_table[talpreds[talrows,]$Var1,"subclass.l1"]
talpreds <- aggregate(talpreds$value,by=list(Var1=talpreds$Var1,Var2=talpreds$Var2),FUN=sum)
talpreds <- as.data.frame(talpreds)
talpreds$cluster <- tals@meta.data[talpreds$Var2,'ctypes_snn_res.0.45']
talpreds$Proportion <- NA
for (cl in unique(talpreds$cluster)){
  talpreds[talpreds$cluster == cl, 'Proportion'] <- talpreds[talpreds$cluster == cl, 'x'] /
    length(unique(talpreds[talpreds$cluster == cl,'Var2']))
}
talpreds <- aggregate(talpreds$Proportion,by=list(Type=talpreds$Var1,Cluster=talpreds$cluster),FUN=sum)
talpreds <- talpreds[talpreds$x > 0,]

lvs <- unique(types_table[!types_table$Subclass.l2 %in% talgroup,'subclass.l1'])
lvs <- c(lvs,talgroup)

talpreds$Type <- factor(talpreds$Type,levels = lvs)


coltable <- read_excel('UCSD_data_032921/Kidney_Reference_Atlas_Summary_Tables.xlsx',
                       sheet='Cluster Color Table',
                       range = 'K1:M101')
coltable <- unique(coltable)
coltable <- coltable[coltable$subclass.l1_label %in% as.character(talpreds$Type),]
coltable <- coltable[coltable$subclass.l1_label != 'FIB',]
cols <- coltable$subclass.l1_color
coltable <- read_excel('UCSD_data_032921/Kidney_Reference_Atlas_Summary_Tables.xlsx',
                       sheet='Cluster Color Table',
                       range = 'H1:J101')
coltable <- unique(coltable)
coltable <- coltable[coltable$subclass.l2_label %in% talgroup,]
cols <- c(cols,coltable$subclass.l2_color)


talothers <- talpreds[talpreds$Type %in% levels(talpreds$Type)[1:13],]

pdf('niches/barplot_all_others.pdf')
ggplot(talothers,aes(x=Cluster,y=x,fill=Type))+
  geom_bar(stat = 'identity',position = 'stack')+
  scale_fill_manual(values = cols)
dev.off()

talpreds$Type <- as.character(talpreds$Type)
talpreds[!talpreds$Type %in% talgroup,'Type'] <- 'Other'
talpreds <- aggregate(talpreds$x,by=list(Type=talpreds$Type,Cluster=talpreds$Cluster),FUN=sum)

talpreds$Type <- factor(talpreds$Type,levels = c('Other',talgroup))

pdf('niches/barplot_tal_sw_relative_cluster_proportions_group_others.pdf')
ggplot(talpreds,aes(x=Cluster,y=x,fill=Type))+
  geom_bar(stat = 'identity',position = 'stack')+
  scale_fill_manual(values = c('#999999',coltable$subclass.l2_color))
dev.off()  




lymph <- levels(talpreds$Type)[c(16:19)]
myel <- levels(talpreds$Type)[c(20:27)]
levs <- c(levels(talpreds$Type)[1:15],'Lymphoid','Myeloid')
talpreds$Type <- as.character(talpreds$Type)
talpreds$Type <- ifelse(talpreds$Type %in% lymph, 'Lymphoid',
                        ifelse(talpreds$Type %in% myel, 'Myeloid', talpreds$Type))
talpreds$Type <- factor(talpreds$Type,levels = levs)
cols <- c('#999999',coltable$subclass.l2_color)
cols <- cols[c(1:15,21,27)]

talpreds <- talpreds[!talpreds$Cluster %in% c(0,4,6,9,10,12),]
talpreds$Cluster <- factor(talpreds$Cluster,levels = c(2,1,8,3,5,7,11))
pdf('niches/barplot_tal_sw_relative_cluster_proportions_group_others_removed_merged.pdf')
ggplot(talpreds,aes(x=Cluster,y=x,fill=Type))+
  geom_bar(stat = 'identity',position = 'stack')+
  scale_fill_manual(values = cols)
dev.off()  



all_merged <- readRDS('all_merged_subclass.l2_norazor.RDS')

all_merged@meta.data$tal_sw_only <- NA
all_merged@meta.data$tal_sw_subclass <- all_merged@meta.data$transfer_subset

all_merged@meta.data[rownames(tals@meta.data),'tal_sw_only'] <- as.character(tals@meta.data$ctypes_snn_res.0.45)
all_merged@meta.data[rownames(tals@meta.data),'tal_sw_subclass'] <- as.character(tals@meta.data$ctypes_snn_res.0.45)


Idents(all_merged) <- all_merged@meta.data$orig.ident
for (sample in c('V19S25-019_XY02-M32','V19S25-019_XY03-M61',
                 'V19S25-016_XY01_18-0006','V19S25-019_XY04-F52',
                 'V10S15-102_XY03_IU-21-015-2','V10S15-102_XY02_IU-21-019-5',
                 'V19S25-017_XY03-13437',#'V19S25-017_XY04-13851',
                 #'V19S25-020_XY03-206850_slide2','V19S25-020_XY04-209400_slide2',
                 'V10S14-085_XY01_20-0038','V10S14-085_XY02_20-0040',
                 'V10S14-085_XY03_21-0056','V10S14-085_XY04_21-0057',
                 'V10S14-086_XY01_21-0055','V10S14-086_XY02_20-0039',
                 'V10S14-086_XY03_21-0063','V10S14-086_XY04_21-0066',
                 'V10S14-087_XY01_21-0061','V10S14-087_XY02_21-0063',
                 'V10S14-087_XY03_21-0064','V10S14-087_XY04_21-0065',
                 'V10S21-388_XY01_21-0068','V10S21-388_XY02_20-0071',
                 'V10S21-388_XY03_20-0072','V10S21-388_XY04_20-0073'
)){
  label <- 'subclass.l2'
  spatial <- readRDS(paste0(label,'/',sample,'_seurat_only.RDS'))
  
  spatial2 <- subset(all_merged,idents = sample)
  
  spots <- unlist(str_split(rownames(spatial2@meta.data),'_'))
  spots <- spots[grep('-',spots)]
  
  spatial@meta.data[spots,'tal_sw_only'] <- spatial2@meta.data$tal_sw_only
  spatial@meta.data[spots,'tal_sw_subclass'] <- spatial2@meta.data$tal_sw_subclass
  maxpred <- spatial@meta.data[,c("tal_sw_only","tal_sw_subclass")]
  maxpred$Barcode <- rownames(spatial@meta.data)
  write.csv(maxpred[,c('Barcode','tal_sw_only')],paste0('niches/',sample,'.tal_sw_only','.csv'),quote=F,row.names=F)
  write.csv(maxpred[,c('Barcode','tal_sw_subclass')],paste0('niches/',sample,'.tal_sw_subclass','.csv'),quote=F,row.names=F)
  
  
}



preds <- tals@assays$predictions[1:74,]
preds <- t(t(preds)/colSums(preds))
preds <- melt(preds)
preds <- as.data.frame(preds)
preds$subclass.l1 <- types_table[as.character(preds$Var1),'subclass.l1']
preds <- preds[preds$subclass.l1 %in% c('TAL','FIB','IMM'),]
preds$talfrac <- NA
for (spot in unique(preds$Var2)){
  preds[preds$Var2 == spot & preds$Var1 == 'aTAL2','talfrac'] <- 
    preds[preds$Var2 == spot & preds$Var1 == 'aTAL2','value'] /
    sum(preds[preds$Var2 == spot & preds$subclass.l1 == 'TAL','value'])
}
preds[preds$subclass.l1 != 'TAL','talfrac'] <- preds[preds$subclass.l1 != 'TAL','value']
preds <- preds[!is.na(preds$talfrac),]
preds[preds$subclass.l1 != 'TAL','subclass.l1'] <- as.character(preds[preds$subclass.l1 != 'TAL','Var1'])
preds <- aggregate(preds$talfrac,by=list(spot=preds$Var2,class=preds$subclass.l1),FUN=sum)
preds <- preds[order(preds$spot),]

preds_table <- dcast(preds,spot ~ class)
ggplot(preds_table,aes(x=`MAC-M2`,y=TAL))+
  geom_point()

talfrac <- preds[preds$class == 'TAL','x']
immfrac <- preds[preds$class == 'IMM','x']
cor(immfrac,talfrac,method = 'pearson')

Idents(all_merged) <- all_merged@meta.data$tal_sw_only
tals_expr <- subset(all_merged,idents = 0:12)
DefaultAssay(tals_expr) <- 'SCT'
saveRDS(tals_expr,'tals_expr.RDS')


markergenes <- read_excel('niches/trajectory_markers.xlsx')
immgenes <- c('COL3A1','LAMA5','SPP1','NRG1','PTPRM',
              'APP','CADM1','EGF','EPHA1','BMP1','NAMPT',
              'SEMA4A','OCLN','CDH1','NEGR1','FGF','SEMA3',
              'THBS','AGRN','GAS','HSPG','NRXN','NECTIN',
              'ALCAM','CNTN','MPZ','KIT','GRN','CXCL')

DotPlot(tals_expr,features = immgenes,col.min = 0)+
  coord_flip()


markers_clean <- c('UMOD','SLC12A1','EGF','WNK1',
                   'REN','CALCR','NOS1','ROBO2','PAPPA2',
                   'SPP1','PROM1','ITGB6','ITGAV','ITGA3',
                   'COL14A1','COL1A2','C7','AEBP1',
                   'IGKC','CCL5','CCL21','FOS','JCHAIN','IGLC2','IGHG3',#'CD38','FCRL5','MS4A1','MZB1',
                   'CD14','CD163','F13A1','HLA-DQA1',
                   'IRAK3','S100A9')
pdf('niches/prelim_markers_short.pdf',height = 4,width = 9)
DotPlot(tals_expr,features = markers_clean,col.min = 0)+
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5))
dev.off()

Idents(tals_expr) <- factor(tals_expr@meta.data$tal_sw_only,levels = c(2,1,8,3,5,7,11,0,4,6,9,10,12))
pdf('niches/neighborhood_markers_short_order.pdf',height = 9.5,width = 3.5)
DotPlot(tals_expr,features = markers_clean,col.min = 0,idents = c(2,1,8,3,5,7,11))+
  scale_y_discrete(labels=c('Medullary (2)','Cortical (1)','Macula Densa (8)','Adaptive (3)',
                           'Stromal (5)','Lymphoid (7)','Myeloid (11)'))+
  ylab('Enrichment of\nTAL Neighborhood')+
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5),
        axis.title.y = element_text(hjust = 0,face = 'bold'))+
  xlab('')
dev.off()

dev.off()


