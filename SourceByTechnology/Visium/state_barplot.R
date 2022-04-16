library(Seurat)
library(ggplot2)
library(reshape2)
library(readxl)


all_merged <- readRDS('all_merged_subclass.l2_norazor.RDS')

coltable <- read_excel('UCSD_data_032921/Kidney_Reference_Atlas_Summary_Tables.xlsx',
                       sheet='Misc Color Table',
                       range = 'Q1:S8')
coltable <- unique(coltable)

types_table <- read_excel('UCSD_data_032921/Kidney_Reference_Atlas_Summary_Tables.xlsx',
                          sheet='Annotations sn10X',range = 'I3:M103')
types_table <- unique(types_table)
rownames(types_table) <- types_table$Subclass.l2

preds <- all_merged@assays$predictions@data[1:74,]
preds <- t(t(preds)/colSums(preds))
preds <- melt(preds)
preds <- as.data.frame(preds)
preds$state <- types_table[as.character(preds$Var1),]$state.l2
preds$condition <- all_merged@meta.data[as.character(preds$Var2),]$condition

bardf <- aggregate(preds$value,by=list(State=preds$state,Condition=preds$condition),FUN=sum)
bardf$Proportion <- NA
for (cond in unique(bardf$Condition)){
  bardf[bardf$Condition == cond,'Proportion'] <- bardf[bardf$Condition == cond,'x']/
    sum(bardf[bardf$Condition == cond,'x'])
}

bardf$State <- factor(bardf$State,levels = coltable$state.l2_label)

coltable <- read_excel('UCSD_data_032921/Kidney_Reference_Atlas_Summary_Tables.xlsx',
                       sheet='Misc Color Table',
                       range = 'H1:J8')
coltable <- unique(coltable)

bardf$Condition <- factor(bardf$Condition,levels = coltable$condition.l1_label)
pdf('newfigs/barplot_proportion_state_condition_by_subclass.pdf',width = 4,height = 4)
ggplot(bardf,aes(x=State,y=Proportion,fill=Condition))+
  geom_bar(stat = 'identity',position = 'dodge')+
  scale_fill_manual(values=coltable$condition.l1_color,
                    labels = c('CKD','AKI','Neph.'))+
  scale_x_discrete(labels=stringr::str_to_title(levels(bardf$State)))
dev.off()  


coltable <- read_excel('UCSD_data_032921/Kidney_Reference_Atlas_Summary_Tables.xlsx',
                       sheet='Misc Color Table',
                       range = 'Q1:S7')
coltable <- unique(coltable)


tst <- matrix(0,nrow = 6,ncol=2,
              dimnames = list(coltable$state.l2_label,c('AKI','CKD')))
for (state in coltable$state.l2_label){
  aki <- c(bardf[bardf$State == state & bardf$Condition == 'AKI','x'],
           sum(bardf[bardf$State != state & bardf$Condition == 'AKI','x']))
  ckd <- c(bardf[bardf$State == state & bardf$Condition == 'CKD','x'],
           sum(bardf[bardf$State != state & bardf$Condition == 'CKD','x']))
  ref <- c(bardf[bardf$State == state & bardf$Condition == 'Ref','x'],
           sum(bardf[bardf$State != state & bardf$Condition == 'Ref','x']))
  print(matrix(c(aki,ref,ncol=2)))
  tst[state,'AKI'] <- fisher.test(matrix(c(aki,ref),ncol=2))$p.value
  tst[state,'CKD'] <- fisher.test(matrix(c(ckd,ref),ncol=2))$p.value
}
write.csv(tst,'newfigs/pvalues_barplot_by_subclass_fisher.csv')

