## prepare cacoa object
library(cacoa)
library(conos)
library(pagoda2)
source("util.func.R")
output.dir <- "./obj"

## functions
#' @param sample.count list of count matrices for each sample
get.sample.count <- function(sample.count, cond=c("AKI", "Ref"), meta=meta.table){
  cond.samples <- unique(meta[meta$condition.l1 %in% cond, ]$specimen)
  return(sample.count[names(sample.count) %in% cond.samples])
}

prepare.cacoa.obj <- function(con, meta=sample.meta, celltype.annot, emb, cond){
  
  sample.group <- setNames(meta[meta$condition.l1 %in% cond, ]$condition.l1, meta[meta$condition.l1 %in% cond, ]$specimen)
  cao <- Cacoa$new(con, sample.groups=sample.group, 
                                  cell.groups = celltype.annot, 
                                  embedding = emb,
                                  n.cores=30, target.level=cond[1], ref.level=cond[2])
  cao$estimateCellLoadings()
  cao$estimateExpressionShiftMagnitudes(min.cells=10, n.cells=1e3, dist="cor", n.subsamples=50)
  #cao$estimateCellDensity()
  cao$estimateCommonExpressionShiftMagnitudes(n.cores=20)
  cao$estimatePerCellTypeDE(independent.filtering=TRUE,  test='DESeq2.Wald',
                            resampling.method='bootstrap', max.resamplings=29)
  cao$estimatePerCellTypeDE(max.cell.count = 50, name='deFixed',
                            resampling.method='bootstrap', max.resamplings=29)
  library(org.Hs.eg.db)
  org <- org.Hs.eg.db
  cao$estimateOntology(type="GO", org.db=org, n.top.genes=500)
  
  
  return(cao)
}


## load data
seurat_tenX <- get(load("./obj/Kidney_KPMP-Biopsy_Ref_10X-R_092020_Seurat_F.rda"))
count.tenX <- seurat_tenX@assays$RNA@counts
meta.table <- seurat_tenX@meta.data
emb <- seurat_tenX@reductions$umap@cell.embeddings

samples <- unique(meta.table$specimen)
sample.count <- lapply(samples, function(r){
  cells <- rownames(meta.table[meta.table$specimen %in% r, ])
  t.count <- count.tenX[, colnames(count.tenX) %in% cells]
  if(ncol(t.count)>200){
    t.count
  }
})
names(sample.count) <- samples
sample.count <- plyr::compact(sample.count)
saveRDS(sample.count, file.path(output.dir, "tenX.sample.count.rds"))

## get AKI vs Ref count
AKI.ref.count <- get.sample.count(sample.count, cond=c("AKI", "Ref"), meta=meta.table)
CKD.ref.count <- get.sample.count(sample.count, cond=c("CKD", "Ref"), meta=meta.table)
AKI.CKD.count <- get.sample.count(sample.count, cond=c("AKI", "CKD"), meta=meta.table)
saveRDS(AKI.ref.count, file.path(output.dir, "tenX.AKI.ref.count.rds"))
saveRDS(CKD.ref.count, file.path(output.dir, "tenX.CKD.ref.count.rds"))
saveRDS(AKI.CKD.count, file.path(output.dir, "tenX.AKI.CKD.count.rds"))

## get conos object
AKI.ref.con <- runConos(AKI.ref.count)
saveRDS(AKI.ref.con, file.path(output.dir, "AKI.ref.con.rds"))

CKD.ref.con <- runConos(CKD.ref.count)
saveRDS(CKD.ref.con, file.path(output.dir, "CKD.ref.con.rds"))

AKI.CKD.con <- runConos(AKI.CKD.count)
saveRDS(AKI.CKD.con, file.path(output.dir, "AKI.CKD.con.rds"))

### prepare cacoa object
sample.meta <- unique(meta.table[, c("specimen", "condition.l1", "subclass.l1")])
celltype.annot <- setNames(meta.table$subclass.l1, rownames(meta.table))

AKI.ref.cao <- prepare.cacoa.obj(AKI.ref.con, meta=sample.meta, celltype.annot, emb, cond=c("AKI", "Ref"))
saveRDS(AKI.ref.cao, file.path(output.dir, "AKI.ref.cao.rds"))

CKD.ref.cao <- prepare.cacoa.obj(CKD.ref.con, meta=sample.meta, celltype.annot, emb, cond=c("CKD", "Ref"))
saveRDS(CKD.ref.cao, file.path(output.dir, "CKD.ref.cao.rds"))

AKI.CKD.cao <- prepare.cacoa.obj(AKI.CKD.con, meta=sample.meta, celltype.annot, emb, cond=c("AKI", "CKD"))
saveRDS(AKI.CKD.cao, file.path(output.dir, "AKI.CKD.cao.rds"))

