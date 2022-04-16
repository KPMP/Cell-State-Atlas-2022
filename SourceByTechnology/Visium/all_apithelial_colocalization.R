library(Seurat)
library(readxl)
library(ggplot2)

all_merged <- readRDS('all_merged_subclass.l2_norazor.RDS')

preds <- all_merged@assays$predictions@data[1:74,]
preds <- t(t(preds)/colSums(preds))
preds_imm <- preds[types_table[types_table$subclass.l1 %in% c('FIB','IMM'),
                               'Subclass.l2'],]
preds_epi <- preds[types_table[types_table$class == 'epithelial cells' ,
                               'Subclass.l2'],]
preds_epi <- melt(preds_epi)
preds_epi <- as.data.frame(preds_epi)
preds_epi$state <- types_table[preds_epi$Var1,'state.l2'] 
preds_epi <- aggregate(preds_epi$value,
                       by=list(spot=preds_epi$Var2,state=preds_epi$state),
                       FUN=sum)
preds_epi <- preds_epi[order(preds_epi$spot),]
preds_epi <- dcast(preds_epi,state ~ spot)
rownames(preds_epi) <- preds_epi$state
preds_epi$state <- NULL
preds_imm <- melt(preds_imm)
preds_imm <- as.data.frame(preds_imm)
preds_imm$Var1 <- as.character(preds_imm$Var1)
preds_imm$imm <- ifelse(preds_imm$Var1 == 'cycMNP','MDC',preds_imm$Var1)
preds_imm <- aggregate(preds_imm$value,
                       by=list(spot=preds_imm$Var2,imm=preds_imm$imm),
                       FUN=sum)
preds_imm <- preds_imm[order(preds_imm$spot),]
preds_imm <- dcast(preds_imm,imm ~ spot)
rownames(preds_imm) <- preds_imm$imm
preds_imm$imm <- NULL
preds <- rbind(preds_epi[c('adaptive - epi','reference','degenerative'),],
               preds_imm)

types_table <- types_table[types_table$Subclass.l2 != 'cycMNP',]
coloc <- matrix(0,nrow = nrow(preds),ncol = ncol(preds),
                dimnames = list(rownames(preds),colnames(preds)))
epilist <- c('adaptive - epi','reference','degenerative')
avgimm <- c()
for (ctype in types_table[types_table$subclass.l1 %in% c('FIB','IMM'),'Subclass.l2']){
  avgimm <- c(avgimm, sum(preds[ctype,preds[ctype,] > 0]) / ncol(preds[,preds[ctype,] > 0]))
}
names(avgimm) <- types_table[types_table$subclass.l1 %in% c('FIB','IMM'),'Subclass.l2']
for(spot in colnames(preds)){
  coloc[max.col(t(preds[epilist,spot])),spot] <- 1
}
for (ctype in types_table[types_table$subclass.l1 %in% c('FIB','IMM'),'Subclass.l2']){
  coloc[ctype,] <- ifelse(t(preds[ctype,])>avgimm[ctype],1,0)
}
write.csv(rowSums(coloc),'niches/coloc_spots_all_epi.csv')

coloc_result <- matrix(0,nrow=length(avgimm),ncol=7,
                       dimnames = list(names(avgimm),c('Coloc. adaptive','adaptive p','adaptive odds',
                                                       'Coloc. degenrative','degenrative p','degenrative odds',
                                                       'Coloc. Reference')))

for (ctype in types_table[types_table$subclass.l1 %in% c('FIB','IMM'),'Subclass.l2']){
  if (length(ncol(coloc[,coloc['adaptive - epi',] == 1 & coloc[ctype,] == 1])) > 0){
    tst <- fisher.test(matrix(c(ncol(coloc[,coloc['adaptive - epi',] == 1 & coloc[ctype,] == 1]),
                                ncol(coloc[,coloc['adaptive - epi',] == 1 & coloc[ctype,] == 0]),
                                ncol(coloc[,coloc['reference',] == 1  & coloc[ctype,] == 1]),
                                ncol(coloc[,coloc['reference',] == 1 & coloc[ctype,] == 0])),ncol = 2))
    coloc_result[ctype,1] <- ncol(coloc[,coloc['adaptive - epi',] == 1 & coloc[ctype,] == 1])
    coloc_result[ctype,2] <- tst$p.value
    coloc_result[ctype,3] <- tst$estimate
  } else{
    coloc_result[ctype,1:3] <- c(0,1,0)
  }
  
  if (length(ncol(coloc[,coloc['degenerative',] == 1 & coloc[ctype,] == 1])) > 0){
    tst <- fisher.test(matrix(c(ncol(coloc[,coloc['degenerative',] == 1 & coloc[ctype,] == 1]),
                                ncol(coloc[,coloc['degenerative',] == 1 & coloc[ctype,] == 0]),
                                ncol(coloc[,coloc['reference',] == 1  & coloc[ctype,] == 1]),
                                ncol(coloc[,coloc['reference',] == 1 & coloc[ctype,] == 0])),ncol = 2))
    coloc_result[ctype,4] <- ncol(coloc[,coloc['degenerative',] == 1 & coloc[ctype,] == 1])
    coloc_result[ctype,5] <- tst$p.value
    coloc_result[ctype,6] <- tst$estimate
  } else{
    coloc_result[ctype,4:6] <- c(0,1,0)
  }
  
  coloc_result[ctype,7] <- ncol(coloc[,coloc['reference',] == 1 & coloc[ctype,] == 1])
}


coltable <- read_excel('UCSD_data_032921/Kidney_Reference_Atlas_Summary_Tables.xlsx',
                       sheet='Misc Color Table',
                       range = 'Q1:S7')
coltable <- unique(coltable)

coloc_plot <- coloc_result[,c(1,4,7)]
colnames(coloc_plot) <- c('Adaptive','Degenerative','Reference')
coloc_plot[,'Adaptive'] <- t(t(coloc_plot[,'Adaptive'])/sum(coloc["adaptive - epi",]))
coloc_plot[,'Degenerative'] <- t(t(coloc_plot[,'Degenerative'])/sum(coloc["degenerative",]))
coloc_plot[,'Reference'] <- t(t(coloc_plot[,'Reference'])/sum(coloc["reference",]))
for (state in rownames(coloc_plot)){
  coloc_plot[state,] <- coloc_plot[state,] / coloc_plot[state,'Reference']
}
coloc_plot <- coloc_plot[!is.na(coloc_plot[,'Reference']),]
coloc_plot <- melt(coloc_plot)

coloc_plot$Var2 <- factor(coloc_plot$Var2,levels = c('Reference','Adaptive','Degenerative'))
coloc_plot$Var1 <- factor(coloc_plot$Var1,levels=c('FIB','aFIB','MYOF','cycMYOF','dFIB','M-FIB','dM-FIB',
                                                   'PL','B','T','NKT','MAST','MAC-M2','MDC','ncMON',
                                                   'cDC','pDC','N'))
pdf('niches/barplot_relative_to_ref_colocalization_events_all_epi.pdf',height = 2.25,width=4.75)
ggplot(coloc_plot,aes(x=Var1,y=value,fill=Var2))+
  geom_bar(stat = 'identity',position = 'dodge')+
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5))+
  xlab('Cell type')+
  ylab('Colocalization')+
  scale_fill_manual(name = "State",values=coltable[c(1,5,2),]$state.l2_color)
dev.off()
