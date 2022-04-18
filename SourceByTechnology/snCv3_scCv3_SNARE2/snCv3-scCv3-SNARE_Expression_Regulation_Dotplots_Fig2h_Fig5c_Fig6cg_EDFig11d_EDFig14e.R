# snCv3/scCv3/SNARE - Expression/Regulation Dotplots (Fig 2h, Fig 5c, Fig 6c&g, ED Fig 11d, ED Fig 14e) ----------------------
library(Signac)
library(Seurat)
library(SeuratDisk)
library(dplyr)
library(Matrix)
library(ggplot2)
library(swne)
library(viridis)
library(chromfunks)

set.seed(1234)
load("color_factors.robj")
load("sc-sn_int_color_factors.robj")

###Prepare seurat objects
load("Kidney_Healthy-Injury_Cell_Atlas_SNARE2-RNA-AC_Seurat_03282022.rda")

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
#Merge to concensus labels 
refquery <- RenameIdents(refquery, 
                         'PT-S1' = 'PT-S1/2',
                         'PT-S2' = 'PT-S1/2',
                         'NKT' = 'NKC/T',
                         'dPT/DTL' = 'dPT')
refquery$subclass.l3 <- Idents(refquery)



###Renal corpuscle dotplot (Fig. 2h)
RQ.GLOM <- subset(refquery, idents = c(
  "POD","PEC","PT-S1/2","MC","EC-GC","EC-AEA","MD","REN","VSMC"))
Idents(RQ.GLOM) <- factor(Idents(RQ.GLOM), levels = c(
  "POD","EC-GC","MC","PEC","PT-S1/2","EC-AEA","VSMC","REN","MD"))

DotPlot(RQ.GLOM, features = c("NPHS2","EMCN","POSTN","VCAM1","SLC5A12","PALMD","MYH11","REN","NOS1"),
        cols = c("light gray", "dodgerblue4"), dot.scale = 8) + RotatedAxis()
rm(RQ.GLOM)


###aSTR dotplot (ED Fig 8f)
RQ.Str <- subset(refquery, idents = c("FIB","dFIB","aFIB","MYOF","cycMYOF"))
Idents(RQ.Str) <- factor(Idents(RQ.Str), levels = c("FIB","dFIB","aFIB","MYOF","cycMYOF"))
ref.st.markers <- c("NEGR1","C7","FBLN5","LAMA2","MEG3","DCN",                                 
                    "B2M","COL1A2","VIM","S100A6","AEBP1","COL1A1","SFRP1","CXCL12","PDGFRA",               #FIB
                    "COL14A1",
                    "POSTN","FLRT2","FGF14","PRRX1","NAV3","IGF1","TWIST2",                    
                    "FAP","SULF1",                                                                           
                    "ACTA2","COL6A3","COL16A1","GLI2","COL5A1","NTM","NOTCH3",
                    "MKI67","TOP2A")

DotPlot(RQ.Str, features = ref.st.markers,
        cols = c("light gray", "dodgerblue4"), dot.scale = 8) + RotatedAxis()
rm(RQ.Str)





###Aging/SASP Gene set scores (Fig. 6g) 
#Prepare objects for score calculation
ref <- KBR
query <- sc.KBR
order <- c(
  "POD","dPOD","PEC","PT-S1/2","PT-S3","aPT","cycPT","dPT",
  "dPT/DTL","DTL2","DTL1","DTL3","dDTL3","ATL","dATL","aTAL1",
  "aTAL2","M-TAL","dM-TAL","C-TAL","dC-TAL","MD","DCT1","DCT2",
  "dDCT","cycDCT","CNT","CNT-PC","dCNT","cycCNT","CCD-PC","OMCD-PC",
  "dOMCD-PC","IMCD","dIMCD","CCD-IC-A","CNT-IC-A","dC-IC-A","OMCD-IC-A",
  "tPC-IC","IC-B","PapE","EC-GC","EC-AEA","EC-DVR","EC-PTC","dEC-PTC",
  "EC-AVR","dEC","cycEC","EC-LYM","MC","REN","VSMC","VSMC/P","dVSMC",
  "MYOF","cycMYOF","FIB","M-FIB","dM-FIB","aFIB","dFIB","B","PL","T",
  "NKC/T","cycNKC/T","MAST","MAC-M2","cycMNP","MDC","cDC","pDC","ncMON","N","SC/NEU")

sen.genes <- read.delim("Scenescence_Marker_Gene_Sets.txt",sep="\t",header=TRUE)

##Tabula.Muris
count.data = GetAssayData(ref, slot = "data")
sen = rownames(count.data)[rownames(count.data) %in% toupper(sen.genes$Tabula.Muris)]
count.data = as.matrix(count.data[rownames(count.data) %in% sen ,])
combined_matrix = rbind(count.data,average=apply(count.data,2,mean))
cor_matrix = cor(t(combined_matrix))
cor_vector = cor_matrix[,dim(cor_matrix)[1]]
count.data_restricted = count.data[rownames(count.data) %in% names(cor_vector[cor_vector >= 0.1]),]
aaa = apply(count.data_restricted,2,mean)
ref$tabula.muris.score <- aaa[rownames(ref@meta.data)]

count.data = GetAssayData(query, slot = "data")
sen = rownames(count.data)[rownames(count.data) %in% toupper(sen.genes$Tabula.Muris)]
count.data = as.matrix(count.data[rownames(count.data) %in% sen ,])
combined_matrix = rbind(count.data,average=apply(count.data,2,mean))
cor_matrix = cor(t(combined_matrix))
cor_vector = cor_matrix[,dim(cor_matrix)[1]]
count.data_restricted = count.data[rownames(count.data) %in% names(cor_vector[cor_vector >= 0.1]),]
aaa = apply(count.data_restricted,2,mean)
query$tabula.muris.score <- aaa[rownames(query@meta.data)]

##Takemon
count.data = GetAssayData(ref, slot = "data")
sen = rownames(count.data)[rownames(count.data) %in% toupper(sen.genes$Takemon)]
count.data = as.matrix(count.data[rownames(count.data) %in% sen ,])
combined_matrix = rbind(count.data,average=apply(count.data,2,mean))
cor_matrix = cor(t(combined_matrix))
cor_vector = cor_matrix[,dim(cor_matrix)[1]]
count.data_restricted = count.data[rownames(count.data) %in% names(cor_vector[cor_vector >= 0.1]),]
aaa = apply(count.data_restricted,2,mean)
ref$Takemon.score <- aaa[rownames(ref@meta.data)]

count.data = GetAssayData(query, slot = "data")
sen = rownames(count.data)[rownames(count.data) %in% toupper(sen.genes$Takemon)]
count.data = as.matrix(count.data[rownames(count.data) %in% sen ,])
combined_matrix = rbind(count.data,average=apply(count.data,2,mean))
cor_matrix = cor(t(combined_matrix))
cor_vector = cor_matrix[,dim(cor_matrix)[1]]
count.data_restricted = count.data[rownames(count.data) %in% names(cor_vector[cor_vector >= 0.1]),]
aaa = apply(count.data_restricted,2,mean)
query$Takemon.score <- aaa[rownames(query@meta.data)]

##Ruscetti
count.data = GetAssayData(ref, slot = "data")
sen = rownames(count.data)[rownames(count.data) %in% toupper(sen.genes$Ruscetti)]
count.data = as.matrix(count.data[rownames(count.data) %in% sen ,])
combined_matrix = rbind(count.data,average=apply(count.data,2,mean))
cor_matrix = cor(t(combined_matrix))
cor_vector = cor_matrix[,dim(cor_matrix)[1]]
count.data_restricted = count.data[rownames(count.data) %in% names(cor_vector[cor_vector >= 0.1]),]
aaa = apply(count.data_restricted,2,mean)
ref$Ruscetti.score <- aaa[rownames(ref@meta.data)]

count.data = GetAssayData(query, slot = "data")
sen = rownames(count.data)[rownames(count.data) %in% toupper(sen.genes$Ruscetti)]
count.data = as.matrix(count.data[rownames(count.data) %in% sen ,])
combined_matrix = rbind(count.data,average=apply(count.data,2,mean))
cor_matrix = cor(t(combined_matrix))
cor_vector = cor_matrix[,dim(cor_matrix)[1]]
count.data_restricted = count.data[rownames(count.data) %in% names(cor_vector[cor_vector >= 0.1]),]
aaa = apply(count.data_restricted,2,mean)
query$Ruscetti.score <- aaa[rownames(query@meta.data)]

##Basisty
count.data = GetAssayData(ref, slot = "data")
sen = rownames(count.data)[rownames(count.data) %in% toupper(sen.genes$Basisty)]
count.data = as.matrix(count.data[rownames(count.data) %in% sen ,])
combined_matrix = rbind(count.data,average=apply(count.data,2,mean))
cor_matrix = cor(t(combined_matrix))
cor_vector = cor_matrix[,dim(cor_matrix)[1]]
count.data_restricted = count.data[rownames(count.data) %in% names(cor_vector[cor_vector >= 0.1]),]
aaa = apply(count.data_restricted,2,mean)
ref$Basisty.score <- aaa[rownames(ref@meta.data)]

count.data = GetAssayData(query, slot = "data")
sen = rownames(count.data)[rownames(count.data) %in% toupper(sen.genes$Basisty)]
count.data = as.matrix(count.data[rownames(count.data) %in% sen ,])
combined_matrix = rbind(count.data,average=apply(count.data,2,mean))
cor_matrix = cor(t(combined_matrix))
cor_vector = cor_matrix[,dim(cor_matrix)[1]]
count.data_restricted = count.data[rownames(count.data) %in% names(cor_vector[cor_vector >= 0.1]),]
aaa = apply(count.data_restricted,2,mean)
query$Basisty.score <- aaa[rownames(query@meta.data)]

##Add to main Refquery object
refquery@meta.data <- refquery@meta.data[c(rownames(ref@meta.data),rownames(query@meta.data)),]
refquery@meta.data$tabula.muris.score <- c(ref@meta.data$tabula.muris.score,query@meta.data$tabula.muris.score)
refquery@meta.data$Takemon.score <- c(ref@meta.data$Takemon.score,query@meta.data$Takemon.score)
refquery@meta.data$Ruscetti.score <- c(ref@meta.data$Ruscetti.score,query@meta.data$Ruscetti.score)
refquery@meta.data$Basisty.score <- c(ref@meta.data$Basisty.score,query@meta.data$Basisty.score)

##Score Plots (ED Fig 14e)
Idents(refquery) <- "subclass.l3"
Idents(refquery) <- factor(Idents(refquery), levels = order)
VlnPlot(refquery, features = c("tabula.muris.score","Takemon.score","Ruscetti.score","Basisty.score"), 
        ncol = 1, pt.size = 0, cols = hsc.l3.cols[levels(Idents(refquery))]) + NoLegend()






###TFBS activities for aPT states (Supplementary Table 20, ED Fig 11d)
mod.assign <- readRDS("pt.snare.module.assignment.rds")
KSAC@meta.data$module <- mod.assign[rownames(KSAC@meta.data)]

Idents(KSAC) <- "module"
KSAC.PT.as <- subset(KSAC, idents = c("black","red","magenta","yellow","turquoise"))
Idents(KSAC.PT.as) <- factor(Idents(KSAC.PT.as), levels = c("black","red","magenta","yellow","turquoise"))

DefaultAssay(KSAC.PT.as) <- "chromvar"
tf.markers <- FindAllMarkers(
  object = KSAC.PT.as,
  max.cells.per.ident = 500,
  only.pos = TRUE,
  test.use = 'LR',
  latent.vars = 'nCount_ATAC'
)

#Subset to TFs enriched and expressed (snCv3/scCv3) 
sn.pt.mod <- readRDS("pt.sn.module.assignment.rds")
sc.pt.mod <- readRDS("pt.sc.module.assignment.rds")
refquery@meta.data[names(sn.pt.mod),]$pt.module <- sn.pt.mod
refquery@meta.data[names(sc.pt.mod),]$pt.module <- sc.pt.mod
Idents(refquery) <- "pt.module"
RQ.PT <- subset(refquery, idents = c("black","red","magenta","yellow","turquoise"))
Idents(RQ.PT) <- factor(Idents(RQ.PT), levels = c("black","red","magenta","yellow","turquoise"))

DefaultAssay(RQ.PT) <- "RNA"
clusters <- levels(Idents(RQ.PT))
rna.counts <- GetAssayData(object = RQ.PT, slot = "counts")
rna.avg.counts <- GetClusterAvg(rna.counts, Idents(RQ.PT), binarize = T)
hist(rna.avg.counts)
quantile(rna.avg.counts)
min.rna.cl.frac <- 0.025
tf.markers$motif <- tf.markers$gene
tf.markers$gene <- gsub("\\s*\\([^\\)]+\\)","",tf.markers$gene) #remove " (var*)"
tf.markers$gene <- toupper(tf.markers$gene)
tf.markers$gene2 <- tf.markers$gene

tf.markers$gene2 <- gsub(":.*","",tf.markers$gene)
tf.markers$gene <- gsub(".*:","",tf.markers$gene)

tf.markers.exp <- lapply(clusters, function(cl) {
  
  cl.expr <- rna.avg.counts[,cl]
  cl.expr <- names(cl.expr[cl.expr > min.rna.cl.frac])
  tf.markers.cl <- tf.markers[tf.markers$cluster == cl,]
  tf.markers.cl <- tf.markers.cl[tf.markers.cl$gene %in% cl.expr | tf.markers.cl$gene2 %in% cl.expr,]
  tf.markers.cl
})
names(tf.markers.exp) <- clusters
tf.markers.exp <- do.call("rbind", lapply(tf.markers.exp, as.data.frame)) 
rownames(tf.markers.exp) <- paste(tf.markers.exp$cluster,tf.markers.exp$motif, sep = ".")
tf.markers.exp

#Plot representative gene expression and TF activities
gene.sub <- c("HNF4A", "THRB", "NR3C2", "MAF", "NFIC", 
              "JUND","JUN","FOSL2","JUNB","FOS","BACH2","ATF2",
              "TGFBR1","TGFBR2","SMAD3","SMAD2",
              "NFKB1","RELB","REL","ATF3","NFKB2","RELA",
              "JAK1","STAT1","STAT2",
              "CREB5","CREB1","CREM",
              "ELF3","ETV6","PAX2","ETS1","EGR1","MEIS1")

DotPlot(RQ.PT, features = gene.sub,
        cols = c("light gray", "dodgerblue4"), dot.scale = 8) + RotatedAxis()

TFBS.sub <- c("HNF4A(var.2)","THRB","NR3C2","MAF","NFIC",
              "JUND(var.2)","JUN","FOSL2","FOS::JUNB","FOS","BACH2","ATF2",
              "Smad2::Smad3",
              "NFKB1","RELB","REL","ATF3","NFKB2","RELA",
              "STAT1::STAT2",
              "Creb5","CREB1","CREM",
              "ELF3","ETV6","Pax2","ETS1","EGR1","MEIS1(var.2)"
)

DefaultAssay(KSAC.PT.as) <- "chromvar"
DotPlot(KSAC.PT.as, features = TFBS.sub, dot.scale = 6, col.min = 0
) + RotatedAxis() + scale_color_viridis() + DarkTheme()







###TFBS activities for aTAL states (Supplementary Table 20, Fig 5c)
snare.mod.assign <- readRDS("tal.snare.module.assignment.rds")
snare.mod.assign <- snare.mod.assign[rownames(KSAC@meta.data)]
KSAC@meta.data$tal.module <- snare.mod.assign

Idents(KSAC) <- "tal.module"
KSAC.TAL <- subset(KSAC, idents = c("black", "pink", "brown", "yellow", "blue"))
Idents(KSAC.TAL) <- factor(Idents(KSAC.TAL), levels = c("black", "pink", "brown", "yellow", "blue"))

DefaultAssay(KSAC.TAL) <- "chromvar"
tf.markers <- FindAllMarkers(
  object = KSAC.TAL,
  max.cells.per.ident = 500,
  only.pos = TRUE,
  test.use = 'LR',
  latent.vars = 'nCount_ATAC'
)

#Subset to TFs enriched and expressed (snCv3/scCv3) 
sn.sc.mod.assign <- readRDS("tal.sn.sc.module.assignment.rds")
sn.sc.mod.assign <- sn.sc.mod.assign[rownames(refquery@meta.data)]

refquery@meta.data$tal.module <- sn.sc.mod.assign
Idents(refquery) <- "tal.module"
RQ.TAL <- subset(refquery, idents = c("black", "pink", "brown", "yellow", "blue"))
Idents(RQ.TAL) <- factor(Idents(RQ.TAL), levels = c("black", "pink", "brown", "yellow", "blue"))

DefaultAssay(RQ.TAL) <- "RNA"
clusters <- levels(Idents(RQ.TAL))
rna.counts <- GetAssayData(object = RQ.TAL, slot = "counts")
rna.avg.counts <- GetClusterAvg(rna.counts, Idents(RQ.TAL), binarize = T)
hist(rna.avg.counts)
quantile(rna.avg.counts)
min.rna.cl.frac <- 0.025
tf.markers$motif <- tf.markers$gene
tf.markers$gene <- gsub("\\s*\\([^\\)]+\\)","",tf.markers$gene) #remove " (var*)"
tf.markers$gene <- toupper(tf.markers$gene)
tf.markers$gene2 <- tf.markers$gene

tf.markers$gene2 <- gsub(":.*","",tf.markers$gene)
tf.markers$gene <- gsub(".*:","",tf.markers$gene)

tf.markers.exp <- lapply(clusters, function(cl) {
  
  cl.expr <- rna.avg.counts[,cl]
  cl.expr <- names(cl.expr[cl.expr > min.rna.cl.frac])
  tf.markers.cl <- tf.markers[tf.markers$cluster == cl,]
  tf.markers.cl <- tf.markers.cl[tf.markers.cl$gene %in% cl.expr | tf.markers.cl$gene2 %in% cl.expr,]
  tf.markers.cl
})
names(tf.markers.exp) <- clusters
tf.markers.exp <- do.call("rbind", lapply(tf.markers.exp, as.data.frame)) 
rownames(tf.markers.exp) <- paste(tf.markers.exp$cluster,tf.markers.exp$motif, sep = ".")
tf.markers.exp


#Plot Representative Subset - RNA and TFBS activities
gene.sub <- c("ESRRG","PPARGC1A","ESRRB","ESRRA",
              "JUND","FOS","JUN","JUNB","BACH1","FOSL2","BACH2","ATF2",
              "ARHGAP29","ROCK2","RAC1","CDC42","RHOA","RHOB","RHOC",
              "STAT3","JAK1","STAT5B","STAT1",
              "SMAD3","SMAD2","TGFBR2","TGFBR1","ITGB6","ITGB8","ITGAV","SKIL",
              "NFKBIZ","NFKB1","RELB","REL",
              "BICC1","DCDC2","TCF7L2","TNIK","PRICKLE1",
              "NR2F2","NR3C1","TFAP2B","KCTD1","ETV6","ZEB1","EHF","ELF3","ETS1","RFX2","CEBPD",
              "FGF13","FGFR1","FGF2",
              "EGF","ERBB4","EGFR","NRG1","ERBB2","PLSCR1","NRG3",
              "MET","MACC1"
              
)
DotPlot(RQ.TAL, features = gene.sub,
        cols = c("light gray", "dodgerblue4"), dot.scale = 8) + RotatedAxis()

TFBS.sub <- c("Esrrg","ESRRB","ESRRA",
              "FOS::JUND","FOSL1::JUN","FOSB::JUNB","BACH1","FOSL2","BACH2","ATF2",
              
              "STAT3","Stat5a::Stat5b","STAT1",
              "Smad2::Smad3",
              
              "REL","RELB","NFKB1",
              "TCF7L2",
              "NR2F2","NR3C1","TFAP2B","ETV6","ZEB1","EHF","ELF3","ETS1","RFX2","CEBPD"
)

DefaultAssay(KSAC.TAL) <- "chromvar"
DotPlot(KSAC.TAL, features = TFBS.sub, dot.scale = 6, col.min = 0
) + RotatedAxis() + scale_color_viridis() + DarkTheme()




###Ligand Receptor Analysis: TAL genes (Fig. 6c)
Idents(refquery) <- "tal.module"
RQ.TAL <- subset(refquery, idents = c("pink", "brown", "yellow", "blue"))
Idents(RQ.TAL) <- factor(Idents(RQ.TAL), levels = c("pink", "brown", "yellow", "blue"))

DotPlot(RQ.TAL, features = c("SLC12A1","EGF","IGF1R","NECTIN2","NRG1","NRG3","BMP6","NRP1","SEMA3C"),
        cols = c("light gray", "dodgerblue4"), dot.scale = 8) + RotatedAxis()




###Ligand Receptor Analysis: STR genes (Fig. 6c)
Idents(refquery) <- "subclass.l3"
RQ.STR <- subset(refquery, idents = c("FIB","aFIB","MYOF"))
Idents(RQ.STR) <- factor(Idents(RQ.STR), levels = c("FIB","aFIB","MYOF"))
DotPlot(RQ.STR, features = c("FBLN5","FLRT2","FGF14","IGF1","IGFBP6","IGFBP5","ERBB4","BMPR2","BMPR1A","ACVR1","BMPR1B",
                             "NRP1","SEMA3C"), 
        cols = c("light gray", "dodgerblue4"), dot.scale = 8) + RotatedAxis()




###Ligand Receptor Analysis: IMM genes (Fig. 6c)
Idents(refquery) <- "subclass.l3"
RQ.IMM <- subset(refquery, idents = c("PL","T","NKC/T","MAC-M2","MDC","ncMON"))
Idents(RQ.IMM) <- factor(Idents(RQ.IMM), levels = c("PL","T","NKC/T","MAC-M2","MDC","ncMON"))
DotPlot(RQ.IMM, features = c("CD226","ERBB4","NRP1"), 
        cols = c("light gray", "dodgerblue4"), dot.scale = 8) + RotatedAxis()



###Ligand Receptor Analysis: TAL genes (Fig. 6c)
RQ.TAL <- subset(RQ.TAL, idents = c("pink", "brown", "yellow", "blue"))
Idents(RQ.TAL) <- factor(Idents(RQ.TAL), levels = c("pink", "brown", "yellow", "blue"))

DotPlot(RQ.TAL, features = c("SLC12A1","EGF","IGF1R","NECTIN2","NRG1","NRG3","BMP6","NRP1","SEMA3C"),
        cols = c("light gray", "dodgerblue4"), dot.scale = 8) + RotatedAxis()




###Aging and SASP expresssion (Fig. 6g)
#PT Modules
sn.pt.mod <- readRDS("pt.sn.module.assignment.rds")
sc.pt.mod <- readRDS("pt.sc.module.assignment.rds")
pt.mod <- c(sn.pt.mod,sc.pt.mod)
cells <- names(pt.mod)
pt.mod <- paste0("PT_", pt.mod)
names(pt.mod) <- cells
refquery@meta.data$module <- NA
refquery@meta.data[names(pt.mod),]$module <- pt.mod

#TAL modules
tal.mod <- readRDS("tal.sn.sc.module.assignment.rds")
cells <- names(tal.mod)
tal.mod <- paste0("TAL_", tal.mod)
names(tal.mod) <- cells
refquery@meta.data[names(tal.mod),]$module <- tal.mod

Idents(refquery) <- "module"
RQ.aEpi <- subset(refquery, idents = c("PT_black","PT_red","PT_magenta","PT_yellow","PT_turquoise",
                                       "TAL_black","TAL_pink","TAL_brown","TAL_yellow","TAL_blue"))
Idents(RQ.aEpi) <- factor(Idents(RQ.aEpi), levels = c("PT_black","PT_red","PT_magenta","PT_yellow","PT_turquoise",
                                                      "TAL_black","TAL_pink","TAL_brown","TAL_yellow","TAL_blue"))
senescent.genes <- c("CDKN1A", "CDKN1B", "CDKN2A")
aEpi <- c("CDH6","HAVCR1","VCAM1","PROM1","ITGB6","CCL2")

DotPlot(RQ.aEpi, features = aEpi,
        cols = c("light gray", "dodgerblue4"), dot.scale = 8) + RotatedAxis()
DotPlot(RQ.aEpi, features = senescent.genes,
        cols = c("light gray", "#99000d"), dot.scale = 8) + RotatedAxis()

VlnPlot(RQ.aEpi, features = c("tabula.muris.score","Takemon.score","Ruscetti.score","Basisty.score"), 
        ncol = 1, pt.size = 0, cols = c("black","red","magenta","yellow","turquoise",
                                        "black","pink","brown","yellow","blue")) + NoLegend()






###Stress TF Activities (Fig. 6g)
#PT Modules
pt.mod <- readRDS("pt.snare.module.assignment.rds")
cells <- names(pt.mod)
pt.mod <- paste0("PT_", pt.mod)
names(pt.mod) <- cells
KSAC@meta.data$module <- NA
KSAC@meta.data[names(pt.mod),]$module <- pt.mod

#TAL modules
tal.mod <- readRDS("tal.snare.module.assignment.rds")
cells <- names(tal.mod)
tal.mod <- paste0("TAL_", tal.mod)
names(tal.mod) <- cells
KSAC@meta.data[names(tal.mod),]$module <- tal.mod

Idents(KSAC) <- "module"
KSAC.as <- subset(KSAC, idents = c("PT_black","PT_red","PT_magenta","PT_yellow","PT_turquoise",
                                   "TAL_black", "TAL_pink", "TAL_brown", "TAL_yellow", "TAL_blue"))
Idents(KSAC.as) <- factor(Idents(KSAC.as), levels = c("PT_black","PT_red","PT_magenta","PT_yellow","PT_turquoise",
                                                      "TAL_black", "TAL_pink", "TAL_brown", "TAL_yellow", "TAL_blue"))
TFBS.sub <- c("XBP1","ATF4","ATF6","HIF1A")

DefaultAssay(KSAC.as) <- "chromvar"
DotPlot(KSAC.as, features = TFBS.sub, dot.scale = 6, col.min = 0
) + RotatedAxis() + scale_color_viridis() + DarkTheme()

##



