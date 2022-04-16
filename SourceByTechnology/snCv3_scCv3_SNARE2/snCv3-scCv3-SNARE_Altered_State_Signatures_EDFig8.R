# snCv3/scCv3/SNARE - Altered State Signatures (ED Fig 8) -------------------------------------------------------------
library(Seurat)
library(SeuratDisk)
library(dplyr)
library(Matrix)
library(ggplot2)

load("sc-sn_int_color_factors.robj")
load("color_factors.robj")

###Prepare seurat objects
KBR <- LoadH5Seurat("Kidney_Healthy-Injury_Cell_Atlas_snCv3_Seurat_03282022.h5Seurat")
Idents(object = KBR) <- "clusters"
Idents(object = KBR) <- factor(Idents(object = KBR), levels = 1:119)
KBR <- subset(KBR, idents = 1:100)
KBR$id <- "ref"
sc.KBR <- LoadH5Seurat("Kidney_Healthy-Injury_Cell_Atlas_scCv3_Seurat_03282022.h5Seurat")
colnames(sc.KBR@meta.data)[colnames(sc.KBR@meta.data) == "state"] <- "state.l1"

#Combined snCv3/scCv3 Object
refquery <- merge(KBR, y = sc.KBR)
umap.coordinates <- rbind(Embeddings(object = KBR, reduction = "umap"),
                          Embeddings(object = sc.KBR, reduction = "ref.umap"))
refquery[["umap"]] <- CreateDimReducObject(embeddings = umap.coordinates, key = "umap_", assay = DefaultAssay(refquery))
refquery <- NormalizeData(refquery, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(refquery)
refquery <- ScaleData(refquery, features = all.genes, split.by = "id", assay = "RNA")
select.markers <- intersect(VariableFeatures(KBR), rownames(refquery))
VariableFeatures(refquery) <- select.markers
Idents(refquery) <- "subclass.l3"


###Plot Figures
###snCv3 - Degenerative State Features (Ext Data Fig 8a)
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



###Matrisome Score Plots (Ext Data Fig 8g)
#core matrisome score
matrisome_set = read.table("https://raw.githubusercontent.com/mahmoudibrahim/KidneyMap/master/assets/public/ecm_genes_human.txt", sep = "\t", header = TRUE)
matrisome = matrisome_set$Gene.Symbol[which(matrisome_set$Division %in% c("Core matrisome"))]
count.data = GetAssayData(KBR, slot = "data")
count.data = as.matrix(count.data[rownames(count.data) %in% toupper(matrisome) ,])
combined_matrix = rbind(count.data,average=apply(count.data,2,mean))
cor_matrix = cor(t(combined_matrix))
cor_vector = cor_matrix[,dim(cor_matrix)[1]]
count.data_restricted = count.data[rownames(count.data) %in% names(cor_vector[cor_vector >= 0.1]),]
aaa = apply(count.data_restricted,2,mean)

KBR$matrisome.score <- aaa[rownames(KBR@meta.data)]

#collagen score
matrisome_set = read.table("https://raw.githubusercontent.com/mahmoudibrahim/KidneyMap/master/assets/public/ecm_genes_human.txt", sep = "\t", header = TRUE)
matrisome = matrisome_set$Gene.Symbol[which(matrisome_set$Category %in% c("Collagens"))]
count.data = GetAssayData(KBR, slot = "data")
count.data = as.matrix(count.data[rownames(count.data) %in% toupper(matrisome) ,])
combined_matrix = rbind(count.data,average=apply(count.data,2,mean))
cor_matrix = cor(t(combined_matrix))
cor_vector = cor_matrix[,dim(cor_matrix)[1]]
count.data_restricted = count.data[rownames(count.data) %in% names(cor_vector[cor_vector >= 0.1]),]
aaa = apply(count.data_restricted,2,mean)

KBR$collagen.score <- aaa[rownames(KBR@meta.data)]

#glycoprotein score
matrisome_set = read.table("https://raw.githubusercontent.com/mahmoudibrahim/KidneyMap/master/assets/public/ecm_genes_human.txt", sep = "\t", header = TRUE)
matrisome = matrisome_set$Gene.Symbol[which(matrisome_set$Category %in% c("ECM Glycoproteins"))]
count.data = GetAssayData(KBR, slot = "data")
count.data = as.matrix(count.data[rownames(count.data) %in% toupper(matrisome) ,])
combined_matrix = rbind(count.data,average=apply(count.data,2,mean))
cor_matrix = cor(t(combined_matrix))
cor_vector = cor_matrix[,dim(cor_matrix)[1]]
count.data_restricted = count.data[rownames(count.data) %in% names(cor_vector[cor_vector >= 0.1]),]
aaa = apply(count.data_restricted,2,mean)

KBR$glycoprotein.score <- aaa[rownames(KBR@meta.data)]

#proteoglycan score
matrisome_set = read.table("https://raw.githubusercontent.com/mahmoudibrahim/KidneyMap/master/assets/public/ecm_genes_human.txt", sep = "\t", header = TRUE)
matrisome = matrisome_set$Gene.Symbol[which(matrisome_set$Category %in% c("Proteoglycans"))]
count.data = GetAssayData(KBR, slot = "data")
count.data = as.matrix(count.data[rownames(count.data) %in% toupper(matrisome) ,])
combined_matrix = rbind(count.data,average=apply(count.data,2,mean))
cor_matrix = cor(t(combined_matrix))
cor_vector = cor_matrix[,dim(cor_matrix)[1]]
count.data_restricted = count.data[rownames(count.data) %in% names(cor_vector[cor_vector >= 0.1]),]
aaa = apply(count.data_restricted,2,mean)

KBR$proteoglycan.score <- aaa[rownames(KBR@meta.data)]

KBR.sub <- subset(KBR, idents = 1:87)
VlnPlot(KBR.sub, features = c("aStr.score", "matrisome.score"), 
        ncol = 1, pt.size = 0, cols = cl.cols[levels(Idents(KBR.sub))]) + NoLegend()





###aEPI Marker Expression (clusters) (ED Fig 8e)

KBR.aEpi <- subset(KBR, idents = c(4:11,26:31,36:37))

Idents(KBR.aEpi) <- factor(Idents(KBR.aEpi), levels = c(4:11,31,36:37,26:27,30,28:29))

mark <- c(   "HAVCR1","VCAM1","PROM1","CST3","CD24",
             "CDHR5","MYO7B","USH1C", 
             "SLC7A8","SLC22A8","SLC22A6","SLC34A1","SLC22A5","SLC36A2","SLC5A12","SLC4A4","SLC13A3",
             "SLC12A6","SLC16A9",	"SLC6A13","SLC22A7","SLC5A10","SLC5A11","SLC6A12",
             "SLC28A1","SLC13A1","SLC16A12","SLC1A1","SLC17A1","SLC2A9","LRP2","CUBN",
             "KCNK13","WNK4",'CLCNKB',"SLC43A2","AKAP6","CPT1A","KCNJ3","SLC25A12",'MPC1',
             "KCNJ1","SLC9A3","SLC16A7","CLCNKA","KCNQ1","CLCN5","ACACB","SLC12A1","UMOD","EGF",
             "MET","LAMC1","ITGAV","NRP1","TNIK","PTPN1","SASH1","COL4A1","ACTN1","ITGA3",
             "DLGAP1","DLG1","CASK","FGFR1","FGFR1OP","FGFR2","FGF13","FGF14",  
             "TGFB1","TGFBR1","TGFBR2","SMAD1","SMAD2","SMAD3","SMAD4",
             "IL15","IL1R1","JAK1","STAT1","STAT3")

DotPlot(KBR.aEpi, features = mark,
        cols = c("light gray", "dodgerblue4"), dot.scale = 8) + RotatedAxis()





###snCv3/scCv3 Stromal Expression plot (Ext Data Fig 8f)
Idents(refquery) <- "subclass.l3"
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





###SNARE2 RNA/AC - Degen State TFBS Activities (ED Fig 8b)
library(Signac)
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(swne)
library(chromfunks)
library(chromVAR)
library(viridis)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)

set.seed(1234)

##Prepare seurat objects
load("Kidney_Healthy-Injury_Cell_Atlas_SNARE2-RNA-AC_Seurat_03282022.rda")

#Combined snCv3/scCv3 Object
refquery <- merge(KBR, y = sc.KBR)
umap.coordinates <- rbind(Embeddings(object = KBR, reduction = "umap"),
                          Embeddings(object = sc.KBR, reduction = "ref.umap"))
refquery[["umap"]] <- CreateDimReducObject(embeddings = umap.coordinates, key = "umap_", assay = DefaultAssay(refquery))
refquery <- NormalizeData(refquery, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(refquery)
refquery <- ScaleData(refquery, features = all.genes, split.by = "id", assay = "RNA")
select.markers <- intersect(VariableFeatures(KBR), rownames(refquery))
VariableFeatures(refquery) <- select.markers


##Degenerative state TF activities
#PTs
Idents(KSAC) <- "subclass.l2"
Idents(KSAC) <- factor(Idents(KSAC), levels = c(
  "POD","dPOD","PEC","PT-S1","PT-S2","PT-S3","aPT","cycPT","dPT","dPT/DTL",
  "DTL2","DTL1","DTL3","dDTL3","ATL","dATL","aTAL1","aTAL2","M-TAL","dM-TAL",
  "C-TAL","dC-TAL","MD","DCT","dDCT","cycDCT","CNT","dCNT","cycCNT","C-PC",
  "M-PC","dM-PC","IMCD","dIMCD","C-IC-A","dC-IC-A","M-IC-A","tPC-IC","IC-B",
  "PapE","EC-GC","EC-AEA","EC-DVR","EC-PTC","dEC-PTC","EC-AVR","dEC","cycEC",
  "EC-LYM","MC","REN","VSMC","VSMC/P","dVSMC","MYOF","cycMYOF","FIB","M-FIB",
  "dM-FIB","aFIB","dFIB","B","PL","T","NKT","MAST","MAC-M2","cycMNP","MDC",
  "cDC","pDC","ncMON","N","SC/NEU"
))
KS.sub.pts13 <- subset(KSAC, idents = c("PT-S1","PT-S2","PT-S3",
                                        "dPT"))
KS.sub.pts13 <- RenameIdents(object = KS.sub.pts13, 
                             "PT-S1" = "PT", 
                             "PT-S2" = "PT",
                             "PT-S3" = "PT")
DefaultAssay(KS.sub.pts13) <- "chromvar"
tf.pts13.markers <- FindAllMarkers(
  object = KS.sub.pts13,
  only.pos = TRUE,
  test.use = 'LR',
  latent.vars = 'nCount_ATAC'
)


#TALs
KS.sub.tal <- subset(KSAC, idents = c("M-TAL","dM-TAL",
                                      "C-TAL","dC-TAL","MD"))
KS.sub.tal <- RenameIdents(object = KS.sub.tal, 
                           "M-TAL" = "TAL", 
                           "dM-TAL" = "dTAL",
                           "C-TAL" = "TAL", 
                           "dC-TAL" = "dTAL",
                           "MD" = "TAL")
DefaultAssay(KS.sub.tal) <- "chromvar"
tf.tal.markers <- FindAllMarkers(
  object = KS.sub.tal,
  only.pos = TRUE,
  test.use = 'LR',
  latent.vars = 'nCount_ATAC'
)


#CNT
KS.sub.cnt <- subset(KSAC, idents = c("CNT","dCNT"))
DefaultAssay(KS.sub.cnt) <- "chromvar"
tf.cnt.markers <- FindAllMarkers(
  object = KS.sub.cnt,
  only.pos = TRUE,
  test.use = 'LR',
  latent.vars = 'nCount_ATAC'
)

#PC
KS.sub.pc <- subset(KSAC, idents = c("C-PC","M-PC","dM-PC","IMCD","dIMCD"))
KS.sub.pc <- RenameIdents(object = KS.sub.pc, 
                          "C-PC" = "PC", 
                          "M-PC" = "PC",
                          "dM-PC" = "dPC", 
                          "IMCD" = "PC",
                          "dIMCD" = "dPC")
DefaultAssay(KS.sub.pc) <- "chromvar"
tf.pc.markers <- FindAllMarkers(
  object = KS.sub.pc,
  only.pos = TRUE,
  test.use = 'LR',
  latent.vars = 'nCount_ATAC'
)

#ICA
KS.sub.ic <- subset(KSAC, idents = c("C-IC-A",  "dC-IC-A",   "M-IC-A"))
KS.sub.ic <- RenameIdents(object = KS.sub.ic, 
                          "C-IC-A" = "IC", 
                          "M-IC-A" = "IC",
                          "dC-IC-A" = "dIC")
DefaultAssay(KS.sub.ic) <- "chromvar"
tf.ic.markers <- FindAllMarkers(
  object = KS.sub.ic,
  only.pos = TRUE,
  test.use = 'LR',
  latent.vars = 'nCount_ATAC'
)

#EC-PTC
KS.sub.ec <- subset(KSAC, idents = c("EC-PTC","dEC-PTC"))
DefaultAssay(KS.sub.ec) <- "chromvar"
tf.ec.markers <- FindAllMarkers(
  object = KS.sub.ec,
  only.pos = TRUE,
  test.use = 'LR',
  latent.vars = 'nCount_ATAC'
)

#FIB
KS.sub.fib <- subset(KSAC, idents = c("FIB","dFIB"))
DefaultAssay(KS.sub.fib) <- "chromvar"
tf.fib.markers <- FindAllMarkers(
  object = KS.sub.fib,
  only.pos = TRUE,
  test.use = 'LR',
  latent.vars = 'nCount_ATAC'
)

#Combine all tables
TF.list = list(tf.pts13.markers,
               tf.tal.markers,
               tf.cnt.markers,
               tf.pc.markers,
               tf.ic.markers,
               tf.ec.markers,
               tf.fib.markers)
celltype = 1:7

cons.TFs <- do.call(rbind, lapply(celltype, function(ct) {
  print(paste("Running for cell type:", ct))
  
  cl.mark <- TF.list[[ct]]
  cl.mark <- cl.mark[cl.mark$p_val < 0.05 & 
                       cl.mark$avg_log2FC > 0.35,]
  return(cl.mark)
  
}))


##Trim TFs to those expressed
Idents(refquery) <- "subclass.l1"
Idents(refquery) <- factor(Idents(refquery), levels = c(
  "POD","PEC","PT","DTL","ATL","TAL","DCT","CNT","PC","IC","PapE",
  "EC","VSM/P","FIB","IMM","NEU"))
refquery <- subset(refquery, idents = c("POD","PT","DTL","ATL","TAL","DCT","CNT","PC","IC",
                                        "EC","VSM/P","FIB"))
Idents(refquery) <- "state.l2"
refquery <- subset(refquery, idents = c("reference","degenerative"))
atac.rna.counts <- GetAssayData(object = refquery, slot = "counts")
atac.rna.idents <- Idents(object = refquery)

sn.cells <- rownames(refquery@meta.data[refquery@meta.data$id == "reference",])
sc.cells <- rownames(refquery@meta.data[refquery@meta.data$id == "query",])

sn.rna.avg.counts <- GetClusterAvg(atac.rna.counts[,sn.cells], atac.rna.idents[sn.cells], binarize = T)
sc.rna.avg.counts <- GetClusterAvg(atac.rna.counts[,sc.cells], atac.rna.idents[sc.cells], binarize = T)

hist(sn.rna.avg.counts)
quantile(sn.rna.avg.counts)
min.rna.cl.frac <- 0.2
sn.degen.expr <- sn.rna.avg.counts[,"degenerative"]
sn.degen.genes <- names(sn.degen.expr[sn.degen.expr > min.rna.cl.frac])
sc.degen.expr <- sc.rna.avg.counts[,"degenerative"]
sc.degen.genes <- names(sc.degen.expr[sc.degen.expr > min.rna.cl.frac])
degen.genes <- unique(c(sn.degen.genes,sc.degen.genes))

cons.TFs$motif <- cons.TFs$gene
cons.TFs$gene <- gsub("\\s*\\([^\\)]+\\)","",cons.TFs$gene)
cons.TFs$gene <- toupper(cons.TFs$gene)
cons.TFs$gene2 <- cons.TFs$gene

cons.TFs$gene2 <- gsub(":.*","",cons.TFs$gene)
cons.TFs$gene <- gsub(".*:","",cons.TFs$gene)

cons.TFs <- cons.TFs[cons.TFs$gene %in% degen.genes | cons.TFs$gene2 %in% degen.genes,]

inj.states <- c("dPT", "dTAL", "dCNT", "dPC", "dIC",
                "dEC-PTC", "dFIB")

inj.TF.counts <- data.frame(cons.TFs[cons.TFs$cluster %in% inj.states,] %>%
                              group_by(motif) %>%
                              tally())
inj.TF.counts <- inj.TF.counts[order(-inj.TF.counts$n),]

inf.TFs <- inj.TF.counts[inj.TF.counts$n >= 3, ]$motif

#Plot ref/degen states
Idents(KSAC) <- "subclass.l1"
Idents(KSAC) <- factor(Idents(KSAC), levels = c(
  "POD","PEC","PT","DTL","ATL","TAL","DCT","CNT","PC","IC","PapE",
  "EC","VSM/P","FIB","IMM","NEU"))
KSAC <- subset(KSAC, idents = c("POD","PT","DTL","ATL","TAL","DCT","CNT","PC","IC",
                                "EC","VSM/P","FIB"))
Idents(KSAC) <- "state.l2"
KSAC <- subset(KSAC, idents = c("reference","degenerative"))
Idents(KSAC) <- "subclass.l2"
Idents(KSAC) <- factor(Idents(KSAC), levels = c(
  "POD","dPOD","PEC","PT-S1","PT-S2","PT-S3","aPT","cycPT","dPT","dPT/DTL",
  "DTL2","DTL1","DTL3","dDTL3","ATL","dATL","aTAL1","aTAL2","M-TAL","dM-TAL",
  "C-TAL","dC-TAL","MD","DCT","dDCT","cycDCT","CNT","dCNT","cycCNT","C-PC",
  "M-PC","dM-PC","IMCD","dIMCD","C-IC-A","dC-IC-A","M-IC-A","tPC-IC","IC-B",
  "PapE","EC-GC","EC-AEA","EC-DVR","EC-PTC","dEC-PTC","EC-AVR","dEC","cycEC",
  "EC-LYM","MC","REN","VSMC","VSMC/P","dVSMC","MYOF","cycMYOF","FIB","M-FIB",
  "dM-FIB","aFIB","dFIB","B","PL","T","NKT","MAST","MAC-M2","cycMNP","MDC",
  "cDC","pDC","ncMON","N","SC/NEU"
))

DefaultAssay(KSAC) <- "chromvar"

DotPlot(KSAC, features = inf.TFs, dot.scale = 6, col.min = 0
) + RotatedAxis() + scale_color_viridis() + DarkTheme()





