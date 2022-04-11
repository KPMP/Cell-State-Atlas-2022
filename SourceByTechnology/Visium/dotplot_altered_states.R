library(Seurat)
library(readxl)
library(ggplot2)
library(reshape2)

all_merged <- readRDS('all_merged_subclass.l2_norazor.RDS')

types_table <- read_excel('UCSD_data_032921/Kidney_Reference_Atlas_Summary_Tables.xlsx',
                          sheet='Annotations sn10X',range = 'I3:M103')
types_table <- as.data.frame(types_table)
types_table <- unique(types_table)
rownames(types_table) <- types_table$Subclass.l2

######
# PT #
######

preds <- as.data.frame(all_merged@assays$predictions@data)[1:74,]
preds <- t(t(preds)/colSums(preds))
preds <- melt(preds)
preds$subclass.l1 <- types_table[as.character(preds$Var1),'subclass.l1']

preds <- aggregate(preds$value,by=list(Spot=preds$Var2,Subclass.l1=preds$subclass.l1),FUN=sum)
preds <- dcast(preds,Subclass.l1 ~ Spot)
rownames(preds) <- preds$Subclass.l1
preds$Subclass.l1 <- NULL
preds <- preds[,rownames(all_merged@meta.data)]
labs <- preds
for (c in 1:ncol(preds)){
  labs[,c] <- rownames(preds)[order(preds[,c],decreasing = T)]
}
all_merged@meta.data$added_subclass <- as.character(labs[1,rownames(all_merged@meta.data)])
Idents(all_merged) <- all_merged@meta.data$added_subclass

preds <- as.data.frame(all_merged@assays$predictions@data)[1:74,]
preds <- t(t(preds)/colSums(preds))
preds <- melt(preds)
preds$subclass.l1 <- paste0(types_table[as.character(preds$Var1),'subclass.l1'],
                            '_',
                            types_table[as.character(preds$Var1),'state.l2'])
preds <- aggregate(preds$value,by=list(Spot=preds$Var2,Subclass.l1=preds$subclass.l1),FUN=sum)
preds <- dcast(preds,Subclass.l1 ~ Spot)
rownames(preds) <- preds$Subclass.l1
preds$Subclass.l1 <- NULL
preds <- preds[,rownames(all_merged@meta.data)]
labs <- preds
for (c in 1:ncol(preds)){
  labs[,c] <- rownames(preds)[order(preds[,c],decreasing = T)]
}

all_merged@meta.data$added_subclass <- as.character(labs[1,rownames(all_merged@meta.data)])
Idents(all_merged) <- all_merged@meta.data$condition
all_merged_con <- subset(all_merged,idents = c('AKI','CKD'))
Idents(all_merged_con) <- all_merged_con@meta.data$added_subclass

pts <- subset(all_merged_con,idents = c("PT_adaptive - epi","PT_cycling",
                                    "PT_degenerative","PT_reference"))
table(pts@meta.data[,c("added_subclass","transfer_subset")])
preds <- as.data.frame(pts@assays$predictions@data)[1:74,]
preds <- t(t(preds)/colSums(preds))
preds <- melt(preds)
preds$subclass.l1 <- paste0(types_table[as.character(preds$Var1),'subclass.l1'],
                            '_',
                            types_table[as.character(preds$Var1),'state.l2'])
sign <- aggregate(preds$value,by=list(Subclass.l2=preds$Var1,Subclass.l1=preds$subclass.l1),FUN=sum)
sign <- sign[grep('PT_',sign$Subclass.l1),]
preds <- aggregate(preds$value,by=list(Spot=preds$Var2,Subclass.l1=preds$subclass.l1),FUN=sum)
preds <- dcast(preds,Subclass.l1 ~ Spot)
rownames(preds) <- preds$Subclass.l1
preds$Subclass.l1 <- NULL
preds <- preds[,rownames(pts@meta.data)]

Idents(pts) <- pts@meta.data$added_subclass
pts <- subset(pts,idents = c("PT_adaptive - epi",
                             "PT_degenerative","PT_reference"))
Idents(pts) <- factor(pts@meta.data$added_subclass,levels=c("PT_adaptive - epi",
                                                            "PT_degenerative","PT_reference"))

final_PT <- c('LRP2','PDZK1','CUBN','GATM','CST3','APOE','ALDOB',
              'SPP1','S100A6','TPM1','CLU','VCAM1','PROM1','ITGB8','HAVCR1')

pdf('dot_alt/alt_PT.pdf',height = 1.9,width = 7)
DotPlot(pts,features = final_PT,
  col.min = 0,col.max = 1,dot.min = 0)+
  theme(axis.text.x = element_text(angle=90,vjust = .5))+
  scale_size_continuous(limits = c(0,100),breaks = c(0,50,100))
dev.off()

#######
# TAL #
#######

tals <- readRDS('tals_expr.RDS')

Idents(tals) <- tals@meta.data$condition
tals <- subset(tals,idents = c('CKD','AKI'))
Idents(tals) <- tals@meta.data$tal_sw_only

preds <- as.data.frame(tals@assays$predictions@data)[1:74,]
preds <- t(t(preds)/colSums(preds))
preds <- preds[types_table[types_table$subclass.l1 == 'TAL',"Subclass.l2"],]
preds <- melt(preds)
preds <- as.data.frame(preds)
preds$state <- types_table[as.character(preds$Var1),"state.l2"]
preds <- aggregate(preds$value,by=list(Spot=preds$Var2,State=preds$state),FUN=sum)
preds <- dcast(preds,State ~ Spot)
rownames(preds) <- preds$State
preds <- preds[,rownames(tals@meta.data)]

labs <- preds
for (c in 1:ncol(preds)){
  labs[,c] <- rownames(preds)[order(preds[,c],decreasing = T)]
}
tals@meta.data$added_subclass <- as.character(labs[1,rownames(tals@meta.data)])
Idents(tals) <- tals@meta.data$added_subclass


adap_TAL <- c("UMOD","SLC12A1","CASR",
              "S100A6", "TPT1", "DEFB1",
              "SPP1","PROM1","ITGB6","CD24")
adap_TAL <- c('UMOD','SLC12A1','DEFB1','S100A6','TPT1','SPP1','CD24','PROM1','ITGB6')
Idents(tals) <- factor(tals@meta.data$added_subclass,levels=c('adaptive - epi','degenerative','reference'))
pdf('dot_alt/alt_TAL_order.pdf',height = 1.9,width = 5.5)
DotPlot(tals,assay='SCT',features = unique(adap_TAL),idents = c('reference','adaptive - epi','degenerative'),
        col.min = 0,col.max = 1,dot.min = 0)+
  theme(axis.text.x = element_text(angle=90,vjust = .5))+
  scale_size_continuous(limits = c(0,100),breaks = c(0,50,100))
dev.off()



#######
# POD #
#######

preds <- as.data.frame(all_merged@assays$predictions@data)[1:74,]
preds <- t(t(preds)/colSums(preds))
preds <- melt(preds)
preds$subclass.l1 <- types_table[as.character(preds$Var1),'subclass.l1']
preds <- aggregate(preds$value,by=list(Spot=preds$Var2,Subclass.l1=preds$subclass.l1),FUN=sum)
preds <- dcast(preds,Subclass.l1 ~ Spot)
rownames(preds) <- preds$Subclass.l1
preds$Subclass.l1 <- NULL

preds <- preds[,rownames(all_merged@meta.data)]
labs <- preds
for (c in 1:ncol(preds)){
  labs[,c] <- rownames(preds)[order(preds[,c],decreasing = T)]
}

all_merged@meta.data$added_subclass <- as.character(labs[1,rownames(all_merged@meta.data)])

all_merged@meta.data$added_subclass <- 'NOT'
all_merged@meta.data[colnames(preds[,preds['POD',]>.2]),"added_subclass"] <- 'keep'

Idents(all_merged) <- all_merged@meta.data$condition
all_merged_con <- subset(all_merged,idents = c('AKI','CKD'))
Idents(all_merged_con) <- all_merged_con@meta.data$added_subclass

pods <- subset(all_merged_con,idents = c('keep'))
table(pods@meta.data[,c("added_subclass","transfer_subset")])
preds <- as.data.frame(pods@assays$predictions@data)[1:74,]
preds <- t(t(preds)/colSums(preds))
preds <- melt(preds)
preds$subclass.l1 <- paste0(types_table[as.character(preds$Var1),'subclass.l1'],
                            '_',
                            types_table[as.character(preds$Var1),'state.l2'])
preds <- aggregate(preds$value,by=list(Spot=preds$Var2,Subclass.l1=preds$subclass.l1),FUN=sum)
preds <- dcast(preds,Subclass.l1 ~ Spot)
rownames(preds) <- preds$Subclass.l1
preds$Subclass.l1 <- NULL
preds <- preds[,rownames(pods@meta.data)]

pod2 <- pods
pod2@meta.data$transfer_subset <- ifelse(pod2@meta.data$transfer_subset == 'EC-GC',
                                        'dPOD',pod2@meta.data$transfer_subset)
Idents(pod2) <- pod2@meta.data$transfer_subset
pod2 <- subset(pod2,idents = c('POD','dPOD'))
Idents(pod2) <- factor(pod2@meta.data$transfer_subset,levels=c('dPOD','POD'))
saveRDS(pod2,'pod2.RDS')
adap_POD <- c('PODXL','NPHS2','CLIC5','SPOCK2','B2M','DC81','S100A6')
pdf('dot_alt/alt_POD_order.pdf',height = 1.7,width = 4.15)
DotPlot(pod2,features = adap_POD,idents = c('POD','dPOD'),
  col.min = 0,col.max = 1,dot.min = 0)+
  theme(axis.text.x = element_text(angle=90,vjust = .5))+
  scale_size_continuous(limits = c(0,100),breaks = c(0,50,100))
dev.off()


#######
# FIB #
#######

preds <- as.data.frame(all_merged@assays$predictions@data)[1:74,]
preds <- t(t(preds)/colSums(preds))
preds <- melt(preds)
preds$subclass.l1 <- paste0(types_table[as.character(preds$Var1),'subclass.l1'],
                            '_',
                            types_table[as.character(preds$Var1),'state.l2'])
preds <- aggregate(preds$value,by=list(Spot=preds$Var2,Subclass.l1=preds$subclass.l1),FUN=sum)
preds <- dcast(preds,Subclass.l1 ~ Spot)
rownames(preds) <- preds$Subclass.l1
preds$Subclass.l1 <- NULL
preds <- preds[,rownames(all_merged@meta.data)]
labs <- preds
for (c in 1:ncol(preds)){
  labs[,c] <- rownames(preds)[order(preds[,c],decreasing = T)]
}


all_merged@meta.data$added_subclass <- as.character(labs[1,rownames(all_merged@meta.data)])
Idents(all_merged) <- all_merged@meta.data$condition
all_merged_con <- subset(all_merged,idents = c('AKI','CKD'))
Idents(all_merged_con) <- all_merged_con@meta.data$added_subclass

fib <- subset(all_merged_con,idents = c("FIB_degenerative","FIB_reference","FIB_adaptive - str"))
table(fib@meta.data[,c("added_subclass","transfer_subset")])
preds <- as.data.frame(fib@assays$predictions@data)[1:74,]
preds <- t(t(preds)/colSums(preds))
preds <- melt(preds)
preds$subclass.l1 <- paste0(types_table[as.character(preds$Var1),'subclass.l1'],
                            '_',
                            types_table[as.character(preds$Var1),'state.l2'])
preds <- aggregate(preds$value,by=list(Spot=preds$Var2,Subclass.l1=preds$subclass.l1),FUN=sum)
preds <- dcast(preds,Subclass.l1 ~ Spot)
rownames(preds) <- preds$Subclass.l1
preds$Subclass.l1 <- NULL
preds <- preds[,rownames(fib@meta.data)]

table(fib@meta.data[,c("added_subclass","transfer_subset")])
table(fib@meta.data[,c("added_subclass","transfer_subset")])
Idents(fib) <- fib@meta.data$added_subclass
Idents(fib) <- factor(fib@meta.data$added_subclass,levels=c("FIB_adaptive - str","FIB_degenerative","FIB_reference"))


final_FIB <- c('FBLN5','PDGFRA','CDH11','IGFBP7','TPT1','S100A6','VIM','ACTA2','COL1A1','COL1A2',
               'KRT19','ABI3BP','PRRX1','FLRT2','IGF1')
pdf('dot_alt/alt_FIB_order.pdf',height = 1.9,width = 7)
DotPlot(fib,features = final_FIB,
  col.min = 0,col.max = 1,dot.min = 0)+
  theme(axis.text.x = element_text(angle=90,vjust = .5))+
  scale_size_continuous(limits = c(0,100),breaks = c(0,50,100))
dev.off()

######
# CD #
######

preds <- as.data.frame(all_merged@assays$predictions@data)[1:74,]
preds <- t(t(preds)/colSums(preds))
preds <- melt(preds)
preds$subclass.l1 <- paste0(types_table[as.character(preds$Var1),'subclass.l1'],
                            '_',
                            types_table[as.character(preds$Var1),'state.l2'])
preds <- aggregate(preds$value,by=list(Spot=preds$Var2,Subclass.l1=preds$subclass.l1),FUN=sum)
preds <- dcast(preds,Subclass.l1 ~ Spot)
rownames(preds) <- preds$Subclass.l1
preds$Subclass.l1 <- NULL
preds <- preds[,rownames(all_merged@meta.data)]
labs <- preds
for (c in 1:ncol(preds)){
  labs[,c] <- rownames(preds)[order(preds[,c],decreasing = T)]
}

all_merged@meta.data$added_subclass <- as.character(labs[1,rownames(all_merged@meta.data)])
Idents(all_merged) <- all_merged@meta.data$condition
all_merged_con <- subset(all_merged,idents = c('AKI','CKD'))
Idents(all_merged_con) <- all_merged_con@meta.data$added_subclass

cds <- subset(all_merged_con,idents = c("PC_degenerative","PC_reference",
                                        "IC_transitioning",
                                        "IC_degenerative","IC_reference"))
table(cds@meta.data[,c("added_subclass","transfer_subset")])
preds <- as.data.frame(cds@assays$predictions@data)[1:74,]
preds <- t(t(preds)/colSums(preds))
preds <- melt(preds)
preds$subclass.l1 <- paste0(types_table[as.character(preds$Var1),'subclass.l1'],
                            '_',
                            types_table[as.character(preds$Var1),'state.l2'])
preds <- aggregate(preds$value,by=list(Spot=preds$Var2,Subclass.l1=preds$subclass.l1),FUN=sum)
preds <- dcast(preds,Subclass.l1 ~ Spot)
rownames(preds) <- preds$Subclass.l1
preds$Subclass.l1 <- NULL
preds <- preds[,rownames(cds@meta.data)]

table(cds@meta.data[,c("added_subclass","transfer_subset")])
table(cds@meta.data[,c("added_subclass","transfer_subset")])
Idents(cds) <- cds@meta.data$added_subclass
Idents(cds) <- factor(cds@meta.data$added_subclass,levels=c("IC_degenerative","PC_degenerative",
                                                            "IC_transitioning",
                                                            "IC_reference","PC_reference"))


final_cds <- c('AQP2','FXYD4','AQP3','GATA3','CALB1','ATP6V0D2',
               'SLC26A7','SLC4A1','SLC4A9','S100A6','S100A11',
               'CLU','MMP7','B2M','FTH1','FTL','APOE')
pdf('dot_alt/alt_CD_order.pdf',height = 2.5,width = 7)
DotPlot(cds,features = final_cds,idents = c("IC_degenerative","PC_degenerative",
                                            "IC_reference","PC_reference"),
  col.min = 0,col.max = 1,dot.min = 0)+
  theme(axis.text.x = element_text(angle=90,vjust = .5))+
  scale_size_continuous(limits = c(0,100),breaks = c(0,50,100))
dev.off()
