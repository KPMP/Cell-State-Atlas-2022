# This script is used to prepare object for Ligand-Receptor analysis

library(CellChat)
library(Seurat)
library(dplyr)
source("./Trajectory/util.func.R")


## load data
sn.seurat <- get(load("./Trajectory/obj/Kidney_KPMP-Biopsy_Ref_10X-R_092020_Seurat_F.rda"))
meta.table <- sn.seurat@meta.data
meta.anno <- meta.table[, c("subclass.l1", "subclass.l2", "state.l1")]
color.table <- read.table("./Trajectory/obj/color.table.txt", sep="\t", comment.char="!", header = T)

tal.cells <- readRDS("./Trajectory/obj/tal.module.assignment.rds") 
meta.anno[rownames(meta.anno) %in% names(tal.cells), ]$subclass.l1 <- tal.cells[match(rownames(meta.anno[rownames(meta.anno) %in% names(tal.cells), ]), names(tal.cells))]
meta.anno <- meta.anno[!is.na(meta.anno$subclass.l1), ]

meta.anno[rownames(meta.anno) %in% names(tal.cells), ]$subclass.l2 <- tal.cells[match(rownames(meta.anno[rownames(meta.anno) %in% names(tal.cells), ]), names(tal.cells))]
meta.anno <- meta.anno[!is.na(meta.anno$subclass.l2), ]

#### TAL to immune cells
### extract tal and immune cells
tal.immune.meta <- meta.anno[meta.anno$subclass.l1 %in% c("IMM", unique(tal.cells)), ]
tal.immune.matrix <- norm.matrix[, colnames(norm.matrix) %in% rownames(tal.immune.meta)]

tal.imm.cellchat <- createCellChat(object=tal.immune.matrix, meta = tal.immune.meta, group.by = "subclass.l3")
tal.imm.cellchat <- addMeta(tal.imm.cellchat, meta = tal.immune.meta)
tal.imm.cellchat <- setIdent(tal.imm.cellchat, ident.use = "subclass.l3")


# using secreted signalling
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
tal.imm.cellchat@DB <- CellChatDB.use
tal.imm.cellchat <- subsetData(tal.imm.cellchat)

tal.imm.cellchat <- tal.imm.cellchat %>% identifyOverExpressedGenes %>% identifyOverExpressedInteractions
tal.imm.cellchat <- projectData(tal.imm.cellchat, PPI.human)

tal.imm.cellchat <- computeCommunProb(tal.imm.cellchat)
tal.imm.cellchat <- filterCommunication(tal.imm.cellchat, min.cells = 10)

# infer signal pathways
tal.imm.cellchat <- computeCommunProbPathway(tal.imm.cellchat) %>% aggregateNet()
saveRDS(tal.imm.cellchat, "./Trajectory/obj/tal.imm.cellchat.obj.rds")