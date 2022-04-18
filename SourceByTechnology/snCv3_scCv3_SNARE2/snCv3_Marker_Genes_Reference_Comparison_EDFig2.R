# snCv3 - Marker Genes and Comparison with Reference Data Sets (ED Fig 2)----------------------------
library(Seurat)
library(plyr)
library(dplyr)
library(Matrix)
library(ggplot2)
library("corrplot")

KBR <- LoadH5Seurat("Kidney_Healthy-Injury_Cell_Atlas_snCv3_Seurat_03282022.h5Seurat")


###Dotplot of marker gene expression (ED Fig 2a)
###Plot a subset of markers from the epi/int plots together
epi.order <- as.character(c(1:3,4:19,20,21,22,24,23,25,26,27,30,31,36,37,39,
                            28,29,32:35,38,40,41,44,45,
                            43,47,42,46,48,49,52,61,56,55,58,60,
                            50,51,53,54,57,59))
int.order <- as.character(c(62:65,67,70,
                            69,66,68,71:74,76:78,80:82,84:85, 
                            79,86,75,83,87,88:99))
KBR.sub <- subset(KBR, idents = c(epi.order,int.order))
Idents(object = KBR.sub) <- factor(Idents(object = KBR.sub), levels = c(epi.order,int.order))


epi.markers <- c("PTPRQ","WT1","NPHS2",                         #POD
                 "CDKN1C", "SPOCK2",                            #dPOD
                 "IGFBP2",                                      #dPOD
                 "CLDN1", "CFH","ALDH1A2",                      #PEC
                 "LRP2",                                        #PT
                 "SLC5A12","SLC22A6",                           #S1/S2
                 "PRODH2","SLC5A2","SLC22A8","SLC7A8",          #S1
                 "SLC34A1","SLC5A10",                           #S2                                  
                 "SLC5A11","SLC7A13",                           #S3
                 
                 "ITGB8","CDH6","HAVCR1","VCAM1",               #aPT
                 "CRYAB","TACSTD2","SLC44A5",                   #TL
                 "AQP1", "UNC5D",                               #DTL1
                 "ADGRL3","ID1",                                #DTL2
                 #DTL3
                 "AKR1B1","SH3GL3",                             #DTL3/ATL
                 "PROX1",                                       #ATL
                 
                 "CASR","SLC12A1","UMOD",                       #TAL
                 "PROM1",                                       #aTAL
                 "ITGB6","FGF13",                               #aTAL
                 "CLDN14", "KCTD16",                            #M-TAL
                 "ANK2",                                        #M-TAL
                 "ESRRB","EGF",                                 #TAL
                 "ENOX1","TMEM207",                             #C-TAL
                 "CLDN16",                                      #C-TAL
                 "NOS1","ROBO2",                                #MD
                 
                 "SLC12A3","TRPM6",                             #DCT
                 #DCT1
                 "SLC8A1","SCN2A","HSD11B2","CALB1",            #CNT
                 "TRPV5",                                       #DCT2
                 "SCNN1G","SCNN1B",                             #CNT                   
                 
                 "GATA3","AQP2","AQP3",                         #PC
                 "PDE10A",                                      #C-PC
                 "KCNK13",                                      #M-PC
                 "FXYD4",                                       #PC
                 "PHACTR1","SLC14A2",                           #IMCD
                 "TP63","GPX2",                                 #PapE
                 "ATP6V0D2",                                    #IC
                 "SLC4A1","SLC26A7",                            #IC-A
                 "KIT","AQP6","CALCA",                          #M-IC-A
                 "SLC4A9","INSRR")                              #IC-B


int.markers <- c("PECAM1","EMCN",                               #EC
                 "HECW2","PLAT","ITGA8",                        #EC-GC
                 "BTNL9","PALMD","AQP1","TM4SF1",               #EC-AEA-DVR
                 "SERPINE2",  #EC-AEA
                 "SLC14A1","ENPP2",                             #EC-DVR
                 "DNASE1L3","CEACAM1",                          #EC-AVR/PTC
                 "PITPNC1","SLCO2A1",                           #EC-PTC
                 "PLVAP","TLL1",                                #EC-AVR
                 "MMRN1","PROX1",                               #EC-LYM
                 
                 "PDGFRB",                                      #VSMC/P
                 
                 "ROBO1",                                       #MC/REN
                 "PIEZO2","POSTN",                              #MC
                 "REN","GRID2",                                 #REN
                 
                 "MYH11","RGS6","MCAM",                         #VSMC               
                 "RGS5","ADGRB3",                               #VSMC/P
                 
                 "COL1A1","COL1A2","C7","DCN",                  #FIB/MyoF
                 
                 "SYNPO2","PCDH7",                              #MyoF
                 "MEG3","LAMA2",                                #FIB/MyoF
                 "COL6A3","GLI2","COL5A1",                      #MyoF
                 "SLC24A3","CDH13","SMOC2","ANO3","RXFP1",      #MyoF
                 
                 "PDGFRA",                                      #FIB
                 
                 "SYT1",                                        #M-FIB
                 "PLCXD3","GREB1L",                             #M-FIB                          
                 "ADAMTSL1","FREM1",                            #M-FIB
                 "MGP","BGN",                                   #M-FIB
                 
                 
                 "ABI3BP","FLRT2", "FGF14","IGF1",              #aFIB
                 
                 "TAGLN","ACTA2",                               #dVSMC/P
                 "PTPRC",                                       #IMM
                 "BANK1","MS4A1",                               #B
                 "IGKC","MZB1",                                 #pB
                 
                 "THEMIS","IL7R","CD96","CD247",                #T
                 "GNLY","NKG7","GZMA",                          #NKT
                 
                 "KIT", "MS4A2",                                #Mast
                 
                 "MRC1","CD163","CD14",                         #MAC
                 "DIAPH3","CENPF","MKI67",                      #cycMAC
                 
                 "MSR1",
                 "ITGAX","HLA-DQA1","CSF2RA",                   #cDC
                 "FLT3","CLEC9A",                               #cDC
                 
                 "IL3RA","CLEC4C",                              #pDC
                 
                 "CTSS","FCN1","FCGR3A",                        #MON
                 "S100A9","S100A8","FCGR3B"                     #NC
                 
)

DotPlot(KBR.sub, features = unique(c(epi.markers,int.markers)),cols = c("lightgrey", "#225ea8"),dot.scale = 4) + RotatedAxis()






###Comparison with reference data sets
###Ransick et al. data (Ex Data Fig 2f) (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE129798)
mcm.meta <- readRDS("/media/Scratch_SSD_Voyager/Blue/KPMP_10X/KPMP_Ref_analysis_12-2020/McMahon_Mouse_2019/pmeta.rds")
mcm.data <- readRDS("/media/Scratch_SSD_Voyager/Blue/KPMP_10X/KPMP_Ref_analysis_12-2020/McMahon_Mouse_2019/expr_data.rds")
mcm.genes <- readRDS("/media/Scratch_SSD_Voyager/Blue/KPMP_10X/KPMP_Ref_analysis_12-2020/McMahon_Mouse_2019/gene_tbl.rds")
mcm.proj <- readRDS("/media/Scratch_SSD_Voyager/Blue/KPMP_10X/KPMP_Ref_analysis_12-2020/McMahon_Mouse_2019/kid_proj.rds")
mcm.otg <- readRDS("/media/Scratch_SSD_Voyager/Blue/KPMP_10X/KPMP_Ref_analysis_12-2020/McMahon_Mouse_2019/otg_tbl.rds")
mcm.otg <- data.frame(mcm.otg)

rownames(mcm.data) <- toupper(rownames(mcm.data))
mcm <- CreateSeuratObject(counts = mcm.data, project = "Mouse Kidney", min.cells = 3, min.features = 200)
mcm.meta <- na.omit(mcm.meta)
colnames(mcm)
mcm <- subset(mcm, cells = rownames(mcm.meta))
mcm <- NormalizeData(mcm, normalization.method = "LogNormalize", scale.factor = 10000)
mcm <- FindVariableFeatures(mcm, selection.method = "vst", nfeatures = 3000)
all.genes <- rownames(mcm)
mcm <- ScaleData(mcm, features = all.genes)

#Use only reference clusters
Idents(KBR) <- "state"
KBR.sub <- subset(KBR, idents = "reference")
Idents(KBR.sub) <- "clusters"
Idents(object = KBR.sub) <- factor(Idents(object = KBR.sub), levels = 1:119)
to.keep <- levels(Idents(KBR.sub))[levels(Idents(KBR.sub)) %in% 1:61]
KBR.sub <- subset(KBR.sub, idents = to.keep)
ave.KBR<-AverageExpression(KBR.sub, features = common.genes, slot = "scale.data")

#reorder Ransick et al. cell types
mcm@meta.data <- mcm.meta[rownames(mcm@meta.data),]
Idents(mcm) <- "Ontology_ID"
order <- c(1:8,"9A","9B",10:20,23,26,28:30,22,25,27,21,24,31,32)
Idents(object = mcm) <- factor(Idents(object = mcm), levels = order)

#Intersect all variable genes
common.genes <- intersect(VariableFeatures(KBR), VariableFeatures(mcm))

ave.KBR<-AverageExpression(KBR.sub, features = common.genes, slot = "scale.data")
ave.mcm<-AverageExpression(mcm, features = common.genes, slot = "scale.data")

ave.cor<-cor(cbind(ave.KBR$RNA,ave.mcm$RNA))
ave.cor<-ave.cor[1:26,27:59]

col2 <- colorRampPalette(c("#053061","#2166AC","#4393C3","#92C5DE","#D1E5F0","#FFFFFF","#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F"))
corrplot(ave.cor, col = col2(500),  method="color",order = "original", cl.lim=c(-1,1))



###Chen et al. scRNA-seq Data (Ex Data Fig 2c) (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE150338)
Idents(object = KBR) <- "clusters"
Idents(object = KBR) <- factor(Idents(object = KBR), levels = 1:119)
KBR.sub <- subset(KBR, idents = 1:100)

load("GSE150338_Seurat_IntegrateData.RData")
Knep.sc <- merge.data.integrated
DefaultAssay(Knep.sc) <- "RNA"
Idents(Knep.sc)
order <- c(10,4,5,13,12,8,0:3,6,9,7,11)
Idents(Knep.sc) <- factor(Idents(Knep.sc), levels = order)
Knep.sc <- ScaleData(Knep.sc)
DefaultAssay(Knep.sc) <- "integrated"

mouse.var <- VariableFeatures(Knep.sc)
DefaultAssay(Knep.sc) <- "RNA"

#Subset to nephron clusters
KBR.sub <- subset(KBR, idents = 26:59)

common.human.genes <- intersect(VariableFeatures(KBR.sub), toupper(mouse.var)) 
common.mouse.genes <- mouse.var[which(toupper(mouse.var) %in% common.human.genes)] 

#Use only reference clusters
ave.KBR <- AverageExpression(KBR.sub, features = common.human.genes, slot = "scale.data")
dim(ave.KBR$RNA)
ave.KBR$RNA <- ave.KBR$RNA[,c("40","41","44","45")]

ave.Kn <- AverageExpression(Knep.sc, features = common.mouse.genes, slot = "scale.data")
rownames(ave.Kn$integrated) <- toupper(rownames(ave.Kn$integrated))
ave.Kn$integrated <- ave.Kn$integrated[common.human.genes,c("0","1","2","3","6","9")]

ave.cor<-cor(cbind(ave.KBR$RNA,ave.Kn$integrated))
ave.cor<-ave.cor[1:4,5:10]
col2 <- colorRampPalette(c("#053061","#2166AC","#4393C3","#92C5DE","#D1E5F0","#FFFFFF","#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F"))
corrplot(ave.cor, col = col2(500),  method="color",order = "original", cl.lim=c(-1,1))




###Chen et al. Bulk Segment Data (Ex Data Fig 2d) (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE150338)
bulk <- read.delim("GSE150338_Mouse_TPM_Replicates.txt",sep="\t",header=TRUE, row.names = 1)
bulk <- bulk[!duplicated(bulk$Gene_symbol),]
rownames(bulk) <- bulk$Gene_symbol
bulk <- bulk[,3:110]

Kn.b <- CreateSeuratObject(counts = bulk, project = "Knepper", min.cells = 3, min.features = 100)
count.data <- GetAssayData(object = Kn.b[["RNA"]], slot = "counts")
Kn.b <- SetAssayData(object = Kn.b,  slot = "data", new.data = count.data, assay = "RNA")
all.genes <- rownames(Kn.b)
Kn.b <- ScaleData(Kn.b, features = all.genes)
Kn.b <- FindVariableFeatures(Kn.b, selection.method = "vst", nfeatures = 3000)
VariableFeatures(Kn.b)

common.human.genes <- intersect(VariableFeatures(KBR.sub), toupper(VariableFeatures(Kn.b))) 
common.mouse.genes <- VariableFeatures(Kn.b)[which(toupper(VariableFeatures(Kn.b)) %in% common.human.genes)] 

#Use only reference clusters
KBR.sub <- subset(KBR, idents = 1:61)

ave.KBR <- AverageExpression(KBR.sub, features = common.human.genes, slot = "scale.data")
order.1 <- c("4","5","6","7","20","21","22","24","31","36",
             "37","39","40","41","44","45","48","49","52","55","56","58","60","61")

ave.Kn <- AverageExpression(Kn.b, features = common.mouse.genes, slot = "scale.data")
rownames(ave.Kn$RNA) <- toupper(rownames(ave.Kn$RNA))
order.2 <- c("PTS1","PTS2","PTS3","DTL1","DTL2","DTL3","ATL","MTAL","CTAL","DCT","CNT","CCD",
             "OMCD","IMCD")

ave.cor<-cor(cbind(ave.KBR$RNA[rownames(ave.Kn$RNA),order.1],ave.Kn$RNA[,order.2]))
ave.cor<-ave.cor[1:24,25:38]

col2 <- colorRampPalette(c("#053061","#2166AC","#4393C3","#92C5DE","#D1E5F0","#FFFFFF","#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F"))
corrplot(ave.cor, col = col2(500),  method="color",order = "original", cl.lim=c(-1,1))



###Kirita et al. data (Ext Data Fig 2e) (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE139107)
meta <- read.table("GSE139107_MouseIRI.metadata.txt",sep="\t",header=TRUE,row.names=1)
H1 <- read.table("GSE139107_MouseIRI_control.dge.txt",sep="\t",header=TRUE,row.names=1)
H2 <- read.table("GSE139107_MouseIRI_4hours.dge.txt",sep="\t",header=TRUE,row.names=1)
H3 <- read.table("GSE139107_MouseIRI_12hours.dge.txt",sep="\t",header=TRUE,row.names=1)
H4 <- read.table("GSE139107_MouseIRI_2days.dge.txt",sep="\t",header=TRUE,row.names=1)
H5 <- read.table("GSE139107_MouseIRI_14days.dge.txt",sep="\t",header=TRUE,row.names=1)
H6 <- read.table("GSE139107_MouseIRI_6weeks.dge.txt",sep="\t",header=TRUE,row.names=1)
H1 <- as.sparse(H1)
H2 <- as.sparse(H2)
H3 <- as.sparse(H3)
H4 <- as.sparse(H4)
H5 <- as.sparse(H5)
H6 <- as.sparse(H6)

merge.sparse = function(listMatrixes) {
  # takes a list of sparse matrixes with different columns and adds them row wise
  
  allRownames <- sort(unique(unlist(lapply(listMatrixes,rownames))))
  allColnames <- as.character(unlist(lapply(listMatrixes,colnames)))
  for (currentMatrix in listMatrixes) {
    newRowLocations <- match(rownames(currentMatrix),allRownames)
    indexes <- which(currentMatrix>0, arr.ind = T)
    newRows <- newRowLocations[indexes[,1]]
    columns <- indexes[,2]
    newMatrix <- sparseMatrix(i=newRows,j=columns, x=currentMatrix@x,
                              dims=c(length(allRownames),max(columns)))
    if (!exists("matrixToReturn")) {
      matrixToReturn <- newMatrix
    }
    else {
      matrixToReturn <- cbind2(matrixToReturn,newMatrix)
    }
  }
  rownames(matrixToReturn) <- allRownames
  colnames(matrixToReturn) <- allColnames
  matrixToReturn  
}
H.list = mget(ls(pattern = "\\H"))
H.data <- merge.sparse(H.list)

rownames(H.data) <- toupper(rownames(H.data))
hki <- CreateSeuratObject(counts = H.data, project = "Mouse Kidney", min.cells = 3, min.features = 200,
                          meta.data = meta)
hki <- NormalizeData(hki, normalization.method = "LogNormalize", scale.factor = 10000)
hki <- FindVariableFeatures(hki, selection.method = "vst", nfeatures = 3000)
all.genes <- rownames(hki)
hki <- ScaleData(hki, features = all.genes)
Idents(hki) <- "celltype"
levels(Idents(hki))
Idents(hki) <- factor(Idents(hki), levels = c("Pod","PEC","PTS1","PTS2","PTS3","NewPT1","NewPT2",
                                              "DTL-ATL","MTAL","CTAL1","CTAL2","MD","DCT","DCT-CNT","CNT",
                                              "PC1","PC2","ICA","ICB","Uro","EC1","EC2","Per","Fib",         
                                              "Tcell","Mø"))
table(Idents(hki))
common.mouse <- intersect(rownames(KBR.sub), VariableFeatures(hki))  
common.human <- intersect(VariableFeatures(KBR.sub)[1:3000], rownames(hki)) 
common.genes <- unique(c(common.mouse, common.human))

#Correlate at the cluster level
Idents(object = KBR.sub) <- "clusters"
Idents(object = KBR.sub) <- factor(Idents(object = KBR.sub), levels = 1:100)

common.mouse <- intersect(rownames(KBR.sub), VariableFeatures(hki))  
common.human <- intersect(VariableFeatures(KBR.sub)[1:3000], rownames(hki)) 
common.genes <- unique(c(common.mouse, common.human))

ave.KBR<-AverageExpression(KBR.sub, features = common.genes, slot = "scale.data")
ave.hki<-AverageExpression(hki, features = common.genes, slot = "scale.data")

ave.cor<-cor(cbind(ave.KBR$RNA,ave.hki$RNA))
ave.cor<-ave.cor[1:12,101:107]

col2 <- colorRampPalette(c("#053061","#2166AC","#4393C3","#92C5DE","#D1E5F0","#FFFFFF","#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F"))
corrplot(ave.cor, col = col2(500),  method="color",order = "original",cl.lim=c(-1,1))




###Kuppe et al. (Ext Data Fig 2g) (https://zenodo.org/record/4059315#.YApbyehKiHs)
Idents(object = KBR) <- "clusters"
Idents(object = KBR) <- factor(Idents(object = KBR), levels = 1:119)
KBR.sub <- subset(KBR, idents = 1:100)

##CD10+ data 
colData <- read.table("Human_CD10plus/PTmap_UMI_counts_colData.txt",sep=",",header=TRUE)
rowData <- read.table("Human_CD10plus/PTmap_UMI_counts_rowData.txt",sep=",",header=TRUE)
Data <- readMM("Human_CD10plus/PTmap_UMI_counts.mtx")
colnames(Data) <- rownames(colData)
rownames(Data) <- rowData$Gene.Symbol

cd10p <- CreateSeuratObject(counts = Data, project = "PT", min.cells = 3, min.features = 200,
                            meta.data = colData)

##CD10- data 
colData <- read.table("CD10negative/kidneyMap_UMI_counts_colData.txt",sep=",",header=TRUE)
rowData <- read.table("CD10negative/kidneyMap_UMI_counts_rowData.txt",sep=",",header=TRUE)
Data <- readMM("CD10negative/kidneyMap_UMI_counts.mtx")
colnames(Data) <- rownames(colData)
rownames(Data) <- rowData$Gene.Symbol

cd10n <- CreateSeuratObject(counts = Data, min.cells = 3, min.features = 200,
                            meta.data = colData)
cd10n$Annotation.Level.3 <- as.character(cd10n$Annotation.Level.3)
cd10n$Annotation.Level.3[cd10n$Annotation.Level.3 == "Fibroblast 6"] <- "Fibroblast 3"
cd10n$Annotation.Level.3[cd10n$Annotation.Level.3 == "Fibroblast 2"] <- "Fibroblast 1"
cd10n$Annotation.Level.3[cd10n$Annotation.Level.3 == "Fibroblast 4"] <- "Fibroblast 2"
table(cd10n$Annotation.Level.3)

#PDGFRb data 
colData <- read.table("Human_PDGFRb/pdgfrbMap_UMI_counts_colData.txt",sep=",",header=TRUE)
rowData <- read.table("Human_PDGFRb/pdgfrbMap_UMI_counts_rowData.txt",sep=",",header=TRUE)
Data <- readMM("Human_PDGFRb/pdgfrbMap_UMI_counts.mtx")
colnames(Data) <- rownames(colData)
rownames(Data) <- rowData$Gene.Symbol

PDGFRb <- CreateSeuratObject(counts = Data, min.cells = 3, min.features = 200,
                             meta.data = colData)
PDGFRb$Annotation.Level.3 <- as.character(PDGFRb$Annotation.Level.3)
PDGFRb$Annotation.Level.3[PDGFRb$Annotation.Level.3 == "Fibroblast 6"] <- "Fibroblast 3"
PDGFRb$Annotation.Level.3[PDGFRb$Annotation.Level.3 == "Fibroblasts 4a"] <- "Fibroblast 2a"
PDGFRb$Annotation.Level.3[PDGFRb$Annotation.Level.3 == "Fibroblast 2"] <- "Fibroblast 1"
PDGFRb$Annotation.Level.3[PDGFRb$Annotation.Level.3 == "Fibroblasts 4b"] <- "Fibroblast 2b"
table(PDGFRb$Annotation.Level.3)


Kr.Seurat <- merge(x = cd10p, y = c(cd10n, PDGFRb))
Idents(Kr.Seurat) <- "Annotation.Level.3"
Kr.Seurat <- NormalizeData(Kr.Seurat, normalization.method = "LogNormalize", scale.factor = 10000)
Kr.Seurat <- FindVariableFeatures(Kr.Seurat, selection.method = "vst", nfeatures = 3000)
all.genes <- rownames(Kr.Seurat)
Kr.Seurat <- ScaleData(Kr.Seurat, features = all.genes)
Idents(Kr.Seurat) <- "Annotation.Level.2"
order <- c("Podocytes","PT","Proximal Tubule","Injured Proximal Tubule",
           "Descending Thin Limb","Thick Ascending Limb","Injured Tubule",
           "Distal Convoluted Tubule","Connecting Tubule","Collecting Duct Principal Cells",
           "Prinicipal Cells","Intercalated Cells","Uroethlial Cells",
           "Glomerular Capillaries","Arterioles","Arterioral Endothelium","Vasa Recta",                   
           "Venular Endothelium","Lymph Endothelium","Injured Endothelium",            
           "Mesangial Cells","Pericytes","Vascular Smooth Muscle Cells",
           "Smooth Muscle Cells","Schwann Cells","Myofibroblast","Myofibroblasts",                   
           "Fibroblast","Fibroblasts",
           "B Cells","Plasma Cells","T Cells","Natural Killer Cells",
           "Macrophages","Dendritic Cells",                                
           "Monocytes","Basophils"                      
)
Idents(Kr.Seurat) <- factor(Idents(Kr.Seurat), levels = order)

#Subset to only the stromal cell types (cortical only)
KBR.st <- subset(KBR.sub, idents = c(71:81,86:87))
Idents(KBR.st) <- factor(Idents(KBR.st), levels = c(71:81,86:87))
Idents(Kr.Seurat) <- "Annotation.Level.3"
levels(Idents(Kr.Seurat))
order <- c("Mesangial Cells 1","Vascular Smooth Muscle Cells","Vascular Smooth Muscle Cells 2",
           "Pericytes","Pericytes 1","Pericytes 2",
           "Myofibroblasts 1","Myofibroblast 1a","Myofibroblast 1b","Myofibroblasts 2a",
           "Myofibroblasts 2b","Myofibroblasts 3a","Myofibroblasts 3b",                                        
           "Fibroblast 1","Fibroblast 2","Fibroblast 2a","Fibroblast 2b","Fibroblast 3")                                                  

Kr.Seurat.st <- subset(Kr.Seurat, idents = order)
Kr.Seurat.st <- FindVariableFeatures(Kr.Seurat.st, selection.method = "vst", nfeatures = 3000)
Idents(Kr.Seurat.st) <- factor(Idents(Kr.Seurat.st), levels = order)

common.genes <- intersect(VariableFeatures(KBR.st), VariableFeatures(Kr.Seurat.st)) 

ave.KBR <- AverageExpression(KBR.st, features = common.genes, slot = "scale.data")
ave.Kr <- AverageExpression(Kr.Seurat.st, features = common.genes, slot = "scale.data")

ave.cor<-cor(cbind(ave.KBR$RNA,ave.Kr$RNA))
ave.cor<-ave.cor[1:13,14:31]

col2 <- colorRampPalette(c("#053061","#2166AC","#4393C3","#92C5DE","#D1E5F0","#FFFFFF","#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F"))
corrplot(ave.cor, col = col2(500),  method="color",order = "original",cl.lim=c(-1,1))




###Immune cell types (Ext Data Fig 2h-i)

#Prepare immune reference data
library(celldex)
library(SingleR)

ref.se <- MonacoImmuneData()

library(scater)
library(Seurat)
library(SeuratDisk)
library(dplyr)
library(Matrix)
library(patchwork)
library(ggplot2)
library("corrplot")

setwd("/media/Scratch_SSD_Voyager/Blue/KPMP_10X/KPMP_Ref_analysis_12-2020")
load("Kidney_KPMP-Biopsy_Ref_10X-R_092020_Seurat_D.rda")
Idents(object = KBR) <- "clusters"
Idents(object = KBR) <- factor(Idents(object = KBR), levels = 1:119)
KBR.sub <- subset(KBR, idents = 1:100)
KBR.imm <- subset(KBR.sub, idents = 88:99)
Idents(KBR.imm) <- factor(Idents(KBR.imm), levels = c(88:99))



##Compare with human sorted immune populations (PBMC) (Ext Data Fig 2i)
logcounts <- assays(ref.se)$logcounts
MID <- CreateSeuratObject(counts = logcounts, min.cells = 0, min.features = 0)
MID@meta.data$label.main <- colData(ref.se)$label.main
MID@meta.data$label.fine <- colData(ref.se)$label.fine
MID@meta.data$label.ont <- colData(ref.se)$label.ont

Idents(MID) <- "label.main"
order <- c("B cells","CD4+ T cells","CD8+ T cells","T cells","NK cells",             
           "Monocytes","Dendritic cells","Neutrophils","Progenitors","Basophils")
Idents(MID) <- factor(Idents(MID), levels = order)

MID <- SetAssayData(
  object = MID,
  slot = "data",
  new.data = logcounts,
  assay = "RNA"
)
all.genes <- rownames(MID)
MID <- ScaleData(MID, features = all.genes)

#Identify and remove poorly aligned subtypes using fine labels
Idents(MID) <- "label.fine"
all.genes <- rownames(MID)

#Use Immune Genes (https://www.immport.org/shared/genelists)
immune.genes <- read.delim("GeneList.txt")
common.genes <- intersect(intersect(rownames(KBR.imm), immune.genes$Symbol),rownames(MID))

#Use only reference clusters
ave.KBR<-AverageExpression(KBR.imm, features = common.genes, slot = "scale.data")
ave.MID<-AverageExpression(MID, features = common.genes, slot = "scale.data")

ave.cor<-cor(cbind(ave.KBR$RNA,ave.MID$RNA))
ave.cor<-ave.cor[1:12,13:41]
ave.cor<-ave.cor[,colSums(ave.cor) > 0.2]

#Remove cell types with poor correlation
MID <- subset(MID, idents = colnames(ave.cor))

Idents(MID) <- "label.main"
order <- c("B cells","CD4+ T cells","CD8+ T cells","T cells","NK cells",             
           "Monocytes","Dendritic cells","Neutrophils","Progenitors","Basophils")
Idents(MID) <- factor(Idents(MID), levels = order)

#Use Immune Genes
common.genes <- intersect(intersect(rownames(KBR.imm), immune.genes$Symbol),rownames(MID))

ave.KBR<-AverageExpression(KBR.imm, features = common.genes, slot = "scale.data")
ave.MID<-AverageExpression(MID, features = common.genes, slot = "scale.data")
ave.cor<-cor(cbind(ave.KBR$RNA,ave.MID$RNA))
ave.cor<-ave.cor[1:12,13:19]
col2 <- colorRampPalette(c("#053061","#2166AC","#4393C3","#92C5DE","#D1E5F0","#FFFFFF","#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F"))
corrplot(ave.cor, col = col2(500),  method="color",order = "original", cl.lim=c(-1,1))



##Compare with mouse immune populations (Ext Data Fig 2h)
library(celldex)
library(SingleR)

ref <- ImmGenData()

logcounts <- assays(ref)$logcounts
rownames(logcounts) <- toupper(rownames(logcounts))

IGD <- CreateSeuratObject(counts = logcounts, min.cells = 0, min.features = 0)
IGD@meta.data$label.main <- colData(ref)$label.main
IGD@meta.data$label.fine <- colData(ref)$label.fine
IGD@meta.data$label.ont <- colData(ref)$label.ont

Idents(IGD) <- "label.main"
order <- c("Epithelial cells","Endothelial cells","Stromal cells","Fibroblasts","Microglia",
           "B cells","B cells, pro","T cells","NKT","NK cells","Mast cells",
           "Macrophages","Monocytes","DC","Neutrophils","Eosinophils","Basophils",             
           "ILC","Tgd","Stem cells")
Idents(IGD) <- factor(Idents(IGD), levels = order)

IGD <- SetAssayData(
  object = IGD,
  slot = "data",
  new.data = logcounts,
  assay = "RNA"
)
all.genes <- rownames(IGD)
IGD <- ScaleData(IGD, features = all.genes)

IGD <- subset(IGD, idents = c("B cells","B cells, pro","T cells","NKT","NK cells","Mast cells",
                              "Macrophages","Monocytes","DC","Neutrophils"))

#Identify and remove poorly aligned subtypes using fine labels
Idents(IGD) <- "label.fine"
all.genes <- rownames(IGD)

#Use Immune Genes (https://www.immport.org/shared/genelists)
immune.genes <- read.delim("GeneList.txt")
common.genes <- intersect(intersect(rownames(KBR.imm), immune.genes$Symbol),rownames(IGD))

ave.KBR<-AverageExpression(KBR.imm, features = common.genes, slot = "scale.data")
ave.IGD<-AverageExpression(IGD, features = common.genes, slot = "scale.data")

ave.cor<-cor(cbind(ave.KBR$RNA,ave.IGD$RNA))
ave.cor<-ave.cor[1:12,13:191]
ave.cor<-ave.cor[,colSums(ave.cor) > 0.2]

#Remove cell types with poor correlation
IGD <- subset(IGD, idents = colnames(ave.cor))

Idents(IGD) <- "label.main"
order <- c("B cells","B cells, pro","T cells","NKT","NK cells","Mast cells",
           "Macrophages","Monocytes","DC","Neutrophils")
Idents(IGD) <- factor(Idents(IGD), levels = order)

common.genes <- intersect(intersect(rownames(KBR.imm), immune.genes$Symbol),rownames(IGD))
ave.KBR<-AverageExpression(KBR.imm, features = common.genes, slot = "scale.data")
ave.IGD<-AverageExpression(IGD, features = common.genes, slot = "scale.data")

ave.cor<-cor(cbind(ave.KBR$RNA,ave.IGD$RNA))
ave.cor<-ave.cor[1:12,13:20]

col2 <- colorRampPalette(c("#053061","#2166AC","#4393C3","#92C5DE","#D1E5F0","#FFFFFF","#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F"))
corrplot(ave.cor, col = col2(500),  method="color",order = "original",
         hclust.method = "ward.D", cl.lim=c(-1,1))




###Predicted Lake et al. 2019 (Nat. Comm.) snDrop-seq Annotations (ED Fig 2b)
##Load snDrop data
#https://drive.google.com/drive/folders/1p-0PNmLDnaTTvz39fr3hA1WKo0rvCWzc?usp=sharing
load("Adult_Normal_Kidney_snDrop_Lake2019_NCOMM_Seuratv3.Robj") 
#meta <- read.table("Kidney_snDrop_Meta.txt",sep="\t",header=TRUE,row.names=1)
#ank.3@meta.data <- meta[rownames(ank.3@meta.data),]
#ank.3$assay <- "snDrop"

##Use variable features as integration anchors
ref.p2 #pagoda2 object generated when clustering snDrop data (not provided)
ref.var.info <- ref.p2$misc$varinfo; ref.var.info <- ref.var.info[order(ref.var.info$lp, decreasing = F),];
ref.var.genes <- rownames(ref.var.info)[1:2000]
query.p2 #pagoda2 object generated when clustering snCv3 data (not provided)
q.var.info <- query.p2$misc$varinfo; q.var.info <- q.var.info[order(q.var.info$lp, decreasing = F),];
query.var.genes <- rownames(q.var.info)[1:2000]

genes.use <- unique(c(ref.var.genes, query.var.genes))
genes.use <- intersect(genes.use, intersect(rownames(ank.3),rownames(KBR))) #2603

NK.anchors <- FindTransferAnchors(reference = ank.3, query = KBR, features = genes.use, 
                                  reduction = "cca")
celltype.predictions <- TransferData(anchorset = NK.anchors, refdata = ank.3$Subclass, 
                                     weight.reduction = "cca")

###

