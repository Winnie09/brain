library(Seurat)
ser <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/humanref/res/integrate.rds')
ref <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/humanref/res/reference.rds')

meta <- ser@meta.data
rownames(meta) <- 
ct <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/integrate/harmony/36nonNormalSeuratGene2000/res/celltype_annotation.rds')  
ser@meta.data <- cbind(meta, celltype=ct[rownames(meta)])

p1 = DimPlot(ref, reduction = "ref.umap", group.by = "Cell.Type", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
p2 = DimPlot(ser, reduction = "ref.umap", group.by = "predicted.celltype.l2", label = TRUE, label.size = 3 ,repel = TRUE) + NoLegend()
p1 + p2

