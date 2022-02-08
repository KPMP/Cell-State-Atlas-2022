library(Seurat)
library(SeuratData)
library(dplyr)
library(data.table)
# read in R object with the post qc-filtered raw read counts from the three TISs
permiereUCSDUCSF = readRDS("Pilot_V1_SC_SN_qcFiltered_rawdata.RDS")

# Seurat integration steps
obj.list <- SplitObject(permiereUCSDUCSF, split.by = "sampletype")

unique(permiereUCSDUCSF@meta.data$sampletype)
obj.list
for (i in 1:length(obj.list)) {
obj.list[[i]] <- NormalizeData(obj.list[[i]], verbose = FALSE)
obj.list[[i]] <- FindVariableFeatures(obj.list[[i]], selection.method = "vst", nfeatures = 2000, verbose = FALSE)
}

tis.anchors <- FindIntegrationAnchors(object.list = obj.list, dims = 1:30, k.anchor = 20)
tis.integrated <- IntegrateData(anchorset = tis.anchors, dims = 1:30)
DefaultAssay(tis.integrated) <- "integrated"
tis.integrated <- ScaleData(tis.integrated, verbose = TRUE )
tis.integrated <- RunPCA(tis.integrated, npcs = 30, verbose = TRUE)
tis.integrated= RunUMAP(tis.integrated,reduction = "pca", dims = 1:18,)
tis.integrated = FindNeighbors(tis.integrated,  dims = 1:18)
tis.integrated= FindClusters(tis.integrated, resolution = .4)
tis.integrated= FindClusters(tis.integrated, resolution = .6)
tis.integrated= FindClusters(tis.integrated, resolution = 1)
DimPlot(tis.integrated, label = T, label.size = 6)

#Find cluster markers
scsn.markers = FindAllMarkers(tis.integrated, only.pos = T)
write.table(scsn.markers, file = "KPMP_PilotV1Integration_02222.txt", sep = '\t', quote = F)


#Assogn cluster names

tis.integrated = RenameIdents(tis.integrated,'0'= 'PT' )
tis.integrated = RenameIdents(tis.integrated,'1'= 'PT' )
tis.integrated = RenameIdents(tis.integrated,'19'= 'PT' )
tis.integrated = RenameIdents(tis.integrated,'24'= 'PT' )
tis.integrated = RenameIdents(tis.integrated,'8'= 'PT' )
tis.integrated = RenameIdents(tis.integrated,'22'= 'PEC' )
tis.integrated = RenameIdents(tis.integrated,'13'= 'DTL' )
tis.integrated = RenameIdents(tis.integrated,'9'= 'ATL/TAL' )
tis.integrated = RenameIdents(tis.integrated,'10'= 'ATL/TAL' )
tis.integrated = RenameIdents(tis.integrated,'28'= 'TAL' )
tis.integrated = RenameIdents(tis.integrated,'6'= 'TAL' )
tis.integrated = RenameIdents(tis.integrated,'7'= 'TAL' )
tis.integrated = RenameIdents(tis.integrated,'12'= 'POD' )
tis.integrated = RenameIdents(tis.integrated,'17'= 'VSMC' )
tis.integrated = RenameIdents(tis.integrated,'14'= 'EC-GC' )
tis.integrated = RenameIdents(tis.integrated,'11'= 'EC-PTC' )
tis.integrated = RenameIdents(tis.integrated,'23'= 'EC-AEA' )
tis.integrated = RenameIdents(tis.integrated,'25'= 'MC' )
tis.integrated = RenameIdents(tis.integrated,'27'= 'FIB' )
tis.integrated = RenameIdents(tis.integrated,'29'= 'B' )
tis.integrated = RenameIdents(tis.integrated,'18'= 'MYL' )
tis.integrated = RenameIdents(tis.integrated,'20'= 'MAC' )
tis.integrated = RenameIdents(tis.integrated,'4'= 'PC' )
tis.integrated = RenameIdents(tis.integrated,'PC'= 'T' )
tis.integrated = RenameIdents(tis.integrated,'21'= 'NKT/NKC' )
tis.integrated = RenameIdents(tis.integrated,'15'= 'DCT' )
tis.integrated = RenameIdents(tis.integrated,'5'= 'CNT' )
tis.integrated = RenameIdents(tis.integrated,'2'= 'PC' )
tis.integrated = RenameIdents(tis.integrated,'26'= 'tPC-IC' )
tis.integrated = RenameIdents(tis.integrated,'16'= 'tPC-IC' )
tis.integrated = RenameIdents(tis.integrated,'3'= 'IC' )


tis.integrated$celltype = Idents(tis.integrated)

save(tis.integrated, file = "scsn.integrated.Robj")

#Plots
my_levels = c('POD','PT','DTL','PEC','ATL/TAL','TAL','DCT','CNT','PC','IC','tPC-IC','EC-AEA','EC-PTC','EC-GC','VSMC','MC','FIB','T','NKT/NKC','MYL','MAC','B')
Idents(tis.integrated) <- factor(Idents(tis.integrated), levels= my_levels)
Idents(tis.integrated) <- 'celltype'
#UMAP
DimPlot(tis.integrated, label = T, label.size = 6, repel = T)
#marker plot
DotPlot(tis.integrated, features = c('NPHS1','NPHS2','GATM','SLC34A1','CUBN','LRP2','CFH','CLDN1', 'VCAM1','SLPI','CLDN3','CLDN10','SLC12A1','EGF','SLC12A3','SLC8A1','STC1','FXYD4','FOXI1','ATP6V0D2','SLC26A7','TAGLN','ACTA2','POSTN','ITGA8','C7','LUM','SERPINE2','EHD3','PLVAP','CD3E','CD3D','NKG7','GZMB','CD14','S100A8','CD68','CD79A','CD79B'),dot.scale = 12)+  theme(text = element_text(size=16)) + theme(axis.text.x = element_text(color = "grey20", size = 16, angle = 90, hjust = .5, vjust = .5, face = "plain"))+ theme(axis.text.y = element_text(color = "grey20", size = 16,  hjust = .5, vjust = .5, face = "plain"))
#nGenes, nCount, mitochondrial percentage plots
VlnPlot(object =tis.integrated, features =c("nFeature_RNA"), ncol =  1 ,pt.size = -1)+  theme(text = element_text(size=22)) + theme(axis.text.x = element_text(color = "grey20", size = 22, angle = 90, hjust = .5, vjust = .5, face = "plain"))+ theme(axis.text.y = element_text(color = "grey20", size = 22,  hjust = .5, vjust = .5, face = "plain"))+theme(legend.position = 'none')
VlnPlot(object =tis.integrated, features =c("nCount_RNA"), ncol =  1 ,pt.size = -1)+  theme(text = element_text(size=22)) + theme(axis.text.x = element_text(color = "grey20", size = 22, angle = 90, hjust = .5, vjust = .5, face = "plain"))+ theme(axis.text.y = element_text(color = "grey20", size = 22,  hjust = .5, vjust = .5, face = "plain"))+theme(legend.position = 'none')
VlnPlot(object =tis.integrated, features =c("percent.mt"), ncol =  1 ,pt.size = -1)+  theme(text = element_text(size=22)) + theme(axis.text.x = element_text(color = "grey20", size = 22, angle = 90, hjust = .5, vjust = .5, face = "plain"))+ theme(axis.text.y = element_text(color = "grey20", size = 22,  hjust = .5, vjust = .5, face = "plain"))+theme(legend.position = 'none')
#ambient mRNA plot
VlnPlot(object = tis.integrated, features =c("LRP2"), ncol =1,  pt.size = 0.00001)
#cell prportion plot
dittoBarPlot(tis.integrated, "celltype", group.by = "sampletype" )+theme(text = element_text(size=20))+ theme(axis.title.x = element_blank(),axis.text.x = element_text(color = "grey20", size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain"))+ theme(axis.text.y = element_text(color = "grey20", size = 20,  hjust = .5, vjust = .5, face = "plain"))+theme( plot.title = element_blank())
dittoBarPlot(tis.integrated, "celltype", group.by = "orig.ident" )+theme(text = element_text(size=20))+ theme(axis.title.x = element_blank(),axis.text.x = element_text(color = "grey20", size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain"))+ theme(axis.text.y = element_text(color = "grey20", size = 20,  hjust = .5, vjust = .5, face = "plain"))+theme( plot.title = element_blank())

