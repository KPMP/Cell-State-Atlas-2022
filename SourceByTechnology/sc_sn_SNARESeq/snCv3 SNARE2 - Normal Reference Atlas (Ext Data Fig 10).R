# snCv3/SNARE2 - Normal Reference Atlas (Ext Data Fig 10) ----------------------------------
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)

int.KCM <- LoadH5Seurat("Kidney_Healthy-Injury_Cell_Atlas_Integrated-snCv3-SNARE-RNA_Seurat_08032021.h5Seurat")
load("color_factors.robj")


###Remove injury state and KPMP disease samples
Idents(int.KCM) <- "state.l1"
int.KCM <- subset(int.KCM, idents = "reference")
Idents(int.KCM) <- "condition.l2"
int.KCM <- subset(int.KCM, idents = "Ref")


###Plots
#Subclass.l3
Idents(object = int.KCM) <- "subclass.l3"
DimPlot(int.KCM, reduction = "umap", label = TRUE, raster=FALSE,
        label.size = 4, repel = TRUE) + ggtitle("Subclass.l3"
        ) + scale_color_manual(values = alpha(sc.l3.cols[levels(Idents(int.KCM))], 0.4), name = "Subclass.l2" 
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2))) + NoLegend()


#Assay
int.KCM$assay[int.KCM$assay == "10X snRNA-Seq"] <- "10X snRNA-seq"
Idents(int.KCM) <- "assay"
assay.cols <- c("#507EB3","gray")
names(assay.cols) <- c("10X snRNA-seq","SNARE-Seq2")

DimPlot(int.KCM, reduction = "umap", label = FALSE, raster=FALSE,
        label.size = 4, repel = TRUE, shuffle = TRUE) + ggtitle("Subclass.l3"
        ) + scale_color_manual(values = alpha(assay.cols[levels(Idents(int.KCM))], 0.4), name = "Subclass.l3" 
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2))) + NoLegend()
assay.cols <- c("gray","#6F3980")
names(assay.cols) <- c("10X snRNA-seq","SNARE-Seq2")

DimPlot(int.KCM, reduction = "umap", label = FALSE, raster=FALSE,
        label.size = 4, repel = TRUE, shuffle = TRUE) + ggtitle("Subclass.l3"
        ) + scale_color_manual(values = alpha(assay.cols[levels(Idents(int.KCM))], 0.4), name = "Subclass.l3" 
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2))) + NoLegend()



#Region.l2
Idents(object = int.KCM) <- "region.l2"
DimPlot(int.KCM, reduction = "umap", label = FALSE, raster=FALSE,
        label.size = 4, repel = TRUE) + ggtitle("Region.l2"
        ) + scale_color_manual(values = alpha(region.l2.cols[levels(Idents(object = int.KCM))], 0.4), name = "Region" 
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2))) 


#Sex
Idents(object = int.KCM) <- "sex"
DimPlot(int.KCM, reduction = "umap", label = FALSE, raster=FALSE,
        label.size = 4, repel = TRUE) + ggtitle("Sex"
        ) + scale_color_manual(values = alpha(sex.cols[levels(Idents(object = int.KCM))], 0.4), name = "Sex" 
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2))) 



###Dotplot of ref markers
Idents(int.KCM) <- "subclass.l3"
Idents(int.KCM) <- factor(Idents(int.KCM), levels = c(
  "POD","dPOD","PEC","PT-S1","PT-S2","PT-S3","aPT","cycPT","dPT","dPT/DTL",
  "DTL2","DTL1","DTL3","dDTL3","ATL","dATL","aTAL1","aTAL2","M-TAL","dM-TAL",
  "C-TAL","dC-TAL","MD","DCT1","DCT2","dDCT","cycDCT","CNT","CNT-PC","dCNT",
  "cycCNT","CCD-PC","OMCD-PC","dOMCD-PC","IMCD","dIMCD","CCD-IC-A","CNT-IC-A",
  "dC-IC-A","OMCD-IC-A","tPC-IC","IC-B","PapE","EC-GC","EC-AEA","EC-DVR","EC-PTC",
  "dEC-PTC","EC-AVR","dEC","cycEC","EC-LYM","MC","REN","VSMC","VSMC/P","dVSMC",
  "MYOF","cycMYOF","FIB","M-FIB","dM-FIB","aFIB","dFIB","B","PL","T","NKT","MAST",
  "MAC-M2","cycMNP","MDC","cDC","pDC","ncMON","N","SC/NEU"
))

markers <- unique(c("PTPRQ",                 #POD
                    "CFH", "ALDH1A2",        #PEC
                    "LRP2",                  #PT
                    "SLC5A12",               #S1/S2
                    "PRODH2",                #S1
                    "SLC34A1",               #S2                                  
                    "SLC5A11","SLC7A13",     #S3
                    "UNC5D","AQP1",          #DTL2
                    "ADGRL3",                #DTL1
                    "SH3GL3",                #DTL3/ATL
                    "PROX1","CLCNKA",        #ATL
                    "UMOD","SLC12A1",        #TAL
                    "ROBO2","NOS1",          #MD
                    "SLC12A3","TRPM6",       #DCT
                    "SLC8A1","HSD11B2",      #CNT
                    "SCNN1G",                #CNT                   
                    "GATA3","AQP2",          #PC
                    "FAM81A",                #M-PC
                    "SLC14A2",               #IMCD
                    "SLC26A7",               #IC-A
                    "LEF1",                  #C-IC-A
                    "KIT","CALCA",           #M-IC-A
                    "SLC4A9","SLC26A4",      #IC-B
                    "TP63","GPX2",           #PAPE
                    "PECAM1",                #EC
                    "HECW2",                 #EC-GC
                    "BTNL9","PALMD",         #EC-AEA-DVR
                    "ENPP2",                 #EC-DVR
                    "DNASE1L3",              #EC-AVR/PTC
                    "PLVAP","GPM6A",         #EC-AVR
                    "MMRN1","CD36",          #EC-LYM
                    "PDGFRB",                #VSMC/P
                    "PIEZO2","POSTN",        #MC
                    "REN",                   #REN
                    "MYH11",                 #C-VSMC/P               
                    "ADGRB3",                #VSMC/P
                    "C7","DCN",              #FIB
                    "PDGFRA",                #FIB
                    "SYT1",                  #M-FIB
                    "BANK1","MS4A1",         #B
                    "MZB1","JCHAIN",         #pB
                    "THEMIS","CD96","CD247", #T
                    "GNLY",                  #NKT
                    "MS4A2",                 #Mast
                    "CD163","CD14",          #MAC
                    "FLT3",                  #cDC
                    "P2RY14","IL3RA",        #pDC
                    "CTSS","FCGR3A",         #MON
                    "S100A9","S100A8","FCGR3B",  #NC
                    "NRXN1","CDH19"          #SC/NEU
))

DotPlot(int.KCM, features = markers,cols = c("lightgrey", "#225ea8"),dot.scale = 4) + RotatedAxis()


