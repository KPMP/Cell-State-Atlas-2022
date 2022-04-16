# snCv3 - Adaptive Epithelial State Gene Set and Score (Supplementary Table 13) -----------------------------------------
library(Seurat)
library(SeuratDisk)
library(dplyr)
library(Matrix)
library(ggplot2)
library(viridis)
library(swne)

KBR <- LoadH5Seurat("Kidney_Healthy-Injury_Cell_Atlas_snCv3_Seurat_03282022.h5Seurat")
Idents(object = KBR) <- "clusters"
Idents(object = KBR) <- factor(Idents(object = KBR), levels = 1:119)
KBR.sub <- subset(KBR, idents = 1:100)


###Conserved aEpi Markers

#Load Function
getconservedDEGs <- function(seurat.obj, celltype) {
  cons.DEGs <- do.call(rbind, lapply(celltype, function(ct) {
    print(paste("Running for cell type:", ct))
    
    mark <- FindConservedMarkers(seurat.obj, ident.1 = ct, grouping.var = "condition.l1",
                                 min.pct = 0.25, only.pos = TRUE, max.cells.per.ident = 300)
    if (length(rownames(mark))>0) {
      mark$gene <- rownames(mark)
      mark$cluster <- rep(ct, length(rownames(mark)))
      rownames(mark) <- paste(mark$cluster, mark$gene, sep = ".")
      return(mark)
    } else {
      print(paste("No DEGs for cell type:", ct))
    }
  }))
  return(cons.DEGs)
}

#Prepare Seurat object
Idents(KBR.sub) <- "subclass.l3"
Idents(KBR.sub) <- factor(Idents(KBR.sub), levels = c(
  "POD","dPOD","PEC","PT-S1","PT-S2","PT-S3","aPT","cycPT","dPT","dPT/DTL",
  "DTL2","DTL1","DTL3","dDTL3","ATL","dATL","aTAL1","aTAL2","M-TAL","dM-TAL",
  "C-TAL","dC-TAL","MD","DCT1","DCT2","dDCT","cycDCT","CNT","CNT-PC","dCNT",
  "cycCNT","CCD-PC","OMCD-PC","dOMCD-PC","IMCD","dIMCD","CCD-IC-A","CNT-IC-A",
  "dC-IC-A","OMCD-IC-A","tPC-IC","IC-B","PapE","EC-GC","EC-AEA","EC-DVR","EC-PTC",
  "dEC-PTC","EC-AVR","dEC","cycEC","EC-LYM","MC","REN","VSMC","VSMC/P","dVSMC",
  "MYOF","cycMYOF","FIB","M-FIB","dM-FIB","aFIB","dFIB","B","PL","T","NKT","MAST",
  "MAC-M2","cycMNP","MDC","cDC","pDC","ncMON","N","SC/NEU"))
table(Idents(KBR.sub))



###Conserved Markers of aPT
KBR.sub.pt <- subset(KBR.sub, idents = c("PT-S1","PT-S2","PT-S3","aPT"))
KBR.sub.pt <- RenameIdents(object = KBR.sub.pt, 
                           "PT-S1" = "PT", 
                           "PT-S2" = "PT",
                           "PT-S3" = "PT")
celltype <- levels(Idents(KBR.sub.pt))
cons.pt.as.DEGs <- getconservedDEGs(KBR.sub.pt, celltype)



###Conserved Markers of aTAL
KBR.sub.tal <- subset(KBR.sub, idents = c("aTAL1","M-TAL","C-TAL"))
KBR.sub.tal <- RenameIdents(object = KBR.sub.tal, 
                            "aTAL1" = "aTAL", 
                            "M-TAL" = "TAL",
                            "C-TAL" = "TAL")
table(Idents(KBR.sub.tal))
Idents(KBR.sub.tal) <- factor(Idents(KBR.sub.tal), levels = c("TAL","aTAL"))
celltype <- levels(Idents(KBR.sub.tal))
cons.tal.as.DEGs <- getconservedDEGs(KBR.sub.tal, celltype)





###Conserved markers across aPT and aTAL
###Combine markers
DEG.list = list(cons.pt.as.DEGs,
                cons.tal.as.DEGs)
celltype = 1:2

cons.DEGs <- do.call(rbind, lapply(celltype, function(ct) {
  print(paste("Running for cell type:", ct))
  
  cl.mark <- DEG.list[[ct]]
  cl.mark <- cl.mark[cl.mark$AKI_p_val < 0.05 & 
                       cl.mark$CKD_p_val < 0.05 & 
                       cl.mark$Ref_p_val < 0.05,]
  top10 <- cl.mark %>% 
    mutate(avg_log2FC = (AKI_avg_log2FC + CKD_avg_log2FC + Ref_avg_log2FC) /3) %>% 
    group_by(cluster) %>% 
    top_n(n = 10, 
          wt = avg_log2FC)
  return(cl.mark)
  
}))




###Select conserved markers that are more specific to the aEpi state 
KBR.sub <- RenameIdents(KBR.sub, "aPT" = "aTub")
KBR.sub <- RenameIdents(KBR.sub, "aTAL1" = "aTub")
KBR.sub <- RenameIdents(KBR.sub, "aTAL2" = "aTub")

aTub.markers <- FindMarkers(KBR.sub, ident.1 = "aTub", features = common.adap.genes,
                            only.pos = TRUE)



###Altered state score
adap.state.genes <- rownames(aTub.markers) #82 genes
count.data = GetAssayData(KBR, slot = "data")
count.data = as.matrix(count.data[rownames(count.data) %in% toupper(adap.state.genes) ,])
combined_matrix = rbind(count.data,average=apply(count.data,2,mean))
cor_matrix = cor(t(combined_matrix))
cor_vector = cor_matrix[,dim(cor_matrix)[1]]
count.data_restricted = count.data[rownames(count.data) %in% names(cor_vector[cor_vector >= 0.1]),]
aaa = apply(count.data_restricted,2,mean)
KBR$aEpi.score <- aaa[rownames(KBR@meta.data)]




###Expanded adaptive -epi state gene set
DEG.list = list(cons.pt.as.DEGs,
                cons.tal.as.DEGs)
celltype = 1:2

cons.DEGs <- do.call(rbind, lapply(celltype, function(ct) {
  print(paste("Running for cell type:", ct))
  
  cl.mark <- DEG.list[[ct]]
  cl.mark <- subset(cl.mark, AKI_avg_log2FC >= 0.6 | CKD_avg_log2FC >= 0.6 | Ref_avg_log2FC >= 0.6)
  cl.mark <- subset(cl.mark, AKI_p_val_adj < 0.05 | CKD_p_val_adj < 0.05 | Ref_p_val_adj < 0.05)
  cl.mark <- cl.mark[cl.mark$max_pval < 0.05,]
  return(cl.mark)
  
}))





###Gene Set Enrichment Analysis (GSEA) (Ext Data Fig. 8c-d)
GSEA = function(gene_list, GO_file, pval) {
  set.seed(54321)
  library(dplyr)
  library(gage)
  library(fgsea)
  library(data.table)
  if ( any( duplicated(names(gene_list)) )  ) {
    warning("Duplicates in gene names")
    gene_list = gene_list[!duplicated(names(gene_list))]
  }
  if  ( !all( order(gene_list, decreasing = TRUE) == 1:length(gene_list)) ){
    warning("Gene list not sorted")
    gene_list = sort(gene_list, decreasing = TRUE)
  }
  myGO = fgsea::gmtPathways(GO_file)
  
  fgRes <- fgsea::fgsea(pathways = myGO, 
                        stats = gene_list,
                        minSize=10,
                        maxSize=500,
                        eps=0) %>% 
    as.data.frame() #%>% 
  #dplyr::filter(padj < !!pval)
  print(dim(fgRes))
  
  ## Filter FGSEA by using gage results. Must be significant and in same direction to keep 
  gaRes = gage::gage(gene_list, gsets=myGO, same.dir=TRUE, set.size =c(10,500))
  
  ups = as.data.frame(gaRes$greater) %>% 
    tibble::rownames_to_column("Pathway") %>% 
    dplyr::filter(!is.na(p.geomean) ) %>%
    dplyr::select("Pathway")
  
  downs = as.data.frame(gaRes$less) %>% 
    tibble::rownames_to_column("Pathway") %>% 
    dplyr::filter(!is.na(p.geomean) ) %>%
    dplyr::select("Pathway")
  
  print(dim(rbind(ups,downs)))
  keepups = as.data.table(fgRes[fgRes$NES > 0 & !is.na(match(fgRes$pathway, ups$Pathway)), ])
  keepdowns = as.data.table(fgRes[fgRes$NES < 0 & !is.na(match(fgRes$pathway, downs$Pathway)), ])
  
  ### Collapse redundant pathways
  Up = fgsea::collapsePathways(keepups[order(pval)], pathways = myGO, stats = gene_list,  nperm = 500, pval.threshold = 0.05)
  Down = fgsea::collapsePathways(keepdowns[order(pval)], myGO, gene_list,  nperm = 500, pval.threshold = 0.05) 
  
  fgRes = fgRes[ !is.na(match(fgRes$pathway, 
                              c( Up$mainPathways, Down$mainPathways))), ] %>% 
    arrange(desc(NES))
  fgRes$pathway = stringr::str_replace(fgRes$pathway, "GO_" , "")
  
  fgRes$Enrichment = ifelse(fgRes$NES > 0, "Up-regulated", "Down-regulated")
  filtRes = rbind(head(fgRes, n = 10),
                  tail(fgRes, n = 10 ))
  g = ggplot(filtRes, aes(reorder(pathway, NES), NES)) +
    geom_segment( aes(reorder(pathway, NES), xend=pathway, y=0, yend=NES)) +
    geom_point( size=5, aes( fill = Enrichment),
                shape=21, stroke=2) +
    scale_fill_manual(values = c("Down-regulated" = "dodgerblue",
                                 "Up-regulated" = "firebrick") ) +
    coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score",
         title="GSEA - Biological Processes") + 
    theme_minimal()
  
  output = list("Results" = fgRes, "Plot" = g)
  return(output)
} #(from https://bioinformaticsbreakdown.com/how-to-gsea/)
GO_file <- "c5.all.v7.2.symbols.gmt"  #GO set (http://www.gsea-msigdb.org/gsea/downloads.jsp#msigdb)


#aPT
S4table <- cons.pt.as.DEGs
S4table$avg_log2FC <- (as.numeric(S4table$AKI_avg_log2FC) + as.numeric(S4table$CKD_avg_log2FC) + as.numeric(S4table$Ref_avg_log2FC)) /3
S4table[S4table$cluster == "aPT",]$avg_log2FC <- as.numeric(S4table[S4table$cluster == "aPT",]$avg_log2FC)
S4table[S4table$cluster == "PT",]$avg_log2FC <- (-1 * as.numeric(S4table[S4table$cluster == "PT",]$avg_log2FC))
gene_list = as.numeric(S4table$avg_log2FC)
names(gene_list) = as.character(S4table$gene)
gene_list = sort(gene_list, decreasing = TRUE)
gene_list = gene_list[!duplicated(names(gene_list))]
res = GSEA(gene_list, GO_file, pval = 0.05)
res$Plot


#aTAL
S4table <- cons.tal.as.DEGs
S4table$avg_log2FC <- (as.numeric(S4table$AKI_avg_log2FC) + as.numeric(S4table$CKD_avg_log2FC) + as.numeric(S4table$Ref_avg_log2FC)) /3
S4table[S4table$cluster == "aTAL",]$avg_log2FC <- as.numeric(S4table[S4table$cluster == "aTAL",]$avg_log2FC)
S4table[S4table$cluster == "TAL",]$avg_log2FC <- (-1 * as.numeric(S4table[S4table$cluster == "TAL",]$avg_log2FC))
gene_list = as.numeric(S4table$avg_log2FC)
names(gene_list) = as.character(S4table$gene)
gene_list = sort(gene_list, decreasing = TRUE)
gene_list = gene_list[!duplicated(names(gene_list))]
res = GSEA(gene_list, GO_file, pval = 0.05)
res$Plot






