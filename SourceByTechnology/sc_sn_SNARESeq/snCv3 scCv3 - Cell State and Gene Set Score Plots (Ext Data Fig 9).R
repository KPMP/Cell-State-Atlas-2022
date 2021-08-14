# snCv3/scCv3 - Cell State and Gene Set Score Plots (Ext Data Fig 9) -------------------------------------------------------------
library(Seurat)
library(SeuratDisk)
library(dplyr)
library(Matrix)
library(ggplot2)

load("sc-sn_int_color_factors.robj")
load("color_factors.robj")

###Prepare seurat objects
refquery <- LoadH5Seurat("Kidney_Healthy-Injury_Cell_Atlas_Integrated-snCv3-scCv3_Seurat_07302021.h5Seurat")
KBR <- LoadH5Seurat("Kidney_Healthy-Injury_Cell_Atlas_snCv3_Seurat_07302021.h5Seurat")
Idents(object = KBR) <- "clusters"
Idents(object = KBR) <- factor(Idents(object = KBR), levels = 1:119)



###snCv3 Cluster Plots (Ext Data Fig 9b)
KBR.sub <- subset(KBR, idents = 1:87)
VlnPlot(KBR.sub, features = c("degen.score","aEpi.score","aStr.score","cyc.score",
                              "matrisome.score"), 
        ncol = 1, pt.size = 0, cols = cl.cols[levels(Idents(KBR.sub))]) + NoLegend()



###scCv3 Subclass Plots (Ext Data Fig 9c)
Idents(object = refquery) <- "id"
sc.KBR <- subset(refquery, idents = "query")
Idents(sc.KBR) <- "subclass.l3"
Idents(sc.KBR) <- factor(Idents(sc.KBR), levels = c(
  "POD","dPOD","PEC","PT-S1/2","PT-S3","aPT","cycPT","dPT",
  "dPT/DTL","DTL2","DTL1","DTL3","dDTL3","ATL","dATL","aTAL1",
  "aTAL2","M-TAL","dM-TAL","C-TAL","dC-TAL","MD","DCT1","DCT2",
  "dDCT","cycDCT","CNT","CNT-PC","dCNT","cycCNT","CCD-PC","OMCD-PC",
  "dOMCD-PC","IMCD","dIMCD","CCD-IC-A","CNT-IC-A","dC-IC-A","OMCD-IC-A",
  "tPC-IC","IC-B","PapE","EC-GC","EC-AEA/DVR","EC-PTC","dEC-PTC",
  "EC-AVR","dEC","cycEC","EC-LYM","MC","REN","VSMC","VSMC/P","dVSMC",
  "MYOF","cycMYOF","FIB","M-FIB","dM-FIB","aFIB","dFIB","B","PL","T",
  "NKC/T","cycNKC/T","MAST","MAC-M2","cycMNP","MDC","cDC","pDC","ncMON","N","SC/NEU"))

sc.KBR <- subset(sc.KBR, idents = c(
  "POD","dPOD","PEC","PT-S1/2","PT-S3","aPT","cycPT","dPT",
  "DTL2","DTL1","DTL3","ATL","aTAL1",
  "aTAL2","M-TAL","dM-TAL","C-TAL","dC-TAL","MD","DCT1","DCT2",
  "dDCT","CNT","CNT-PC","dCNT","CCD-PC","OMCD-PC",
  "dOMCD-PC","IMCD","CCD-IC-A","CNT-IC-A","dC-IC-A","OMCD-IC-A",
  "tPC-IC","IC-B","PapE","EC-GC","EC-AEA/DVR","EC-PTC",
  "EC-AVR","dEC","cycEC","EC-LYM","MC","VSMC","VSMC/P","dVSMC",
  "MYOF","FIB","M-FIB","aFIB","dFIB"))

VlnPlot(sc.KBR, features = c("degen.score","aEpi.score","aStr.score","cyc.score",
                             "matrisome.score"), 
        ncol = 1, pt.size = 0, cols = hsc.l3.cols[levels(Idents(sc.KBR))]) + NoLegend()




###snCv3 - Degenerative State Features (Ext Data Fig 9d)
KBR.sub <- subset(KBR, idents = 1:100)
Idents(object = KBR.sub) <- "state.l2"
KBR.sub.inj <- subset(KBR.sub, idents = c("reference","degenerative"))
Idents(object = KBR.sub.inj) <- "state.l1"
Idents(object = KBR.sub.inj) <- factor(Idents(object = KBR.sub.inj), levels = c("reference","altered"))
KBR.sub.inj$state.l1 <- Idents(object = KBR.sub.inj)

Idents(KBR.sub.inj) <- "subclass.l1"
Idents(object = KBR.sub.inj) <- factor(Idents(object = KBR.sub.inj), 
                                       levels = c("POD","PEC","PT","DTL","ATL","TAL",
                                                  "DCT","CNT","PC","IC","PapE","EC","VSM/P",
                                                  "FIB","IMM","NEU"))
KBR.sub.inj <- subset(KBR.sub.inj, idents = c("POD","PT","DTL","ATL","TAL",
                                              "DCT","CNT","PC","IC","EC","VSM/P",
                                              "FIB"))


VlnPlot(object = KBR.sub.inj, features = c("IGFBP7","CLU","CST3","degen.score","nFeature_RNA","percent.er","percent.mt"), 
        ncol = 1, pt.size = -1, split.by = "state.l1", cols = as.character(state.l1.cols))



###scCv3 - Degenerative State Features (Ext Data Fig 9e)
Idents(object = refquery) <- "id"
sc.KBR <- subset(refquery, idents = "query")
Idents(object = sc.KBR) <- "state.l2"
sc.KBR.inj <- subset(sc.KBR, idents = c("reference","degenerative"))
Idents(object = sc.KBR.inj) <- "state.l1"
Idents(object = sc.KBR.inj) <- factor(Idents(object = sc.KBR.inj), levels = c("reference","altered"))
sc.KBR.inj$state.l1 <- Idents(object = sc.KBR.inj)

Idents(sc.KBR.inj) <- "subclass.l1"
Idents(object = sc.KBR.inj) <- factor(Idents(object = sc.KBR.inj), 
                                      levels = c("POD","PEC","PT","DTL","ATL","TAL",
                                                 "DCT","CNT","PC","IC","PapE","EC","VSM/P",
                                                 "FIB","IMM","NEU"))
sc.KBR.inj <- subset(sc.KBR.inj, idents = c("POD","PT","TAL",
                                            "DCT","CNT","PC","IC","EC","VSM/P",
                                            "FIB"))

VlnPlot(object = sc.KBR.inj, features = c("IGFBP7","CLU","CST3","degen.score","nFeature_RNA","percent.er","percent.mt"), 
        ncol = 1, pt.size = -1, split.by = "state.l1", cols = as.character(state.l1.cols)) 



###snCv3/scCv3 Score plots split by condition - Ext Data Fig 9g
refquery <- LoadH5Seurat("Kidney_Healthy-Injury_Cell_Atlas_Integrated-snCv3-scCv3_Seurat_07302021.h5Seurat")
VlnPlot(refquery, features = c("aEpi.score","cyc.score"), split.by = "condition.l1",
        ncol = 1, pt.size = 0, cols = as.character(cond.l1.cols[c("AKI","CKD","Ref")])) + NoLegend()



###snCv3/scCv3 Stromal Expression plot (Ext Data Fig 9a)
RQ.Str <- subset(refquery, idents = c("FIB","dFIB","aFIB","MYOF","cycMYOF"))
Idents(RQ.Str) <- factor(Idents(RQ.Str), levels = c("FIB","dFIB","aFIB","MYOF","cycMYOF"))
ref.st.markers <- c("NEGR1","C7","FBLN5","LAMA2","MEG3","DCN",
                    "B2M","COL1A2","VIM","S100A6","AEBP1","COL1A1","SFRP1","CXCL12","PDGFRA",
                    "COL14A1",
                    "POSTN","FLRT2","FGF14","PRRX1","NAV3","IGF1","TWIST2",
                    "FAP","SULF1",
                    "ACTA2","COL6A3","COL16A1","GLI2","COL5A1","NTM","NOTCH3",
                    "MKI67","TOP2A")

DotPlot(RQ.Str, features = ref.st.markers,
        cols = c("light gray", "dodgerblue4"), dot.scale = 8) + RotatedAxis()


