library(Seurat)
ser <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/humanref/res/integrate.rds')
ref <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/humanref/res/reference.rds')

meta <- ser@meta.data

ct <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/integrate/harmony/36nonNormalSeuratGene2000/res/celltype_annotation.rds')  
ser@meta.data <- cbind(meta, celltype=ct[rownames(meta)])

p1 = DimPlot(ref, reduction = "umap", group.by = "Cell.Type", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
p2 = DimPlot(ser, reduction = "ref.umap", group.by = "celltype", label = TRUE, label.size = 3 ,repel = TRUE) + NoLegend()
p1 + p2

png('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/humanref/plot/alltogether_cortex_v3_celltype2.png', width = 1800, height = 1400, res = 250)
DimPlot(reference, reduction = "umap", group.by = "Cell.Type", label = F, label.size = 3 ,repel = TRUE) #+ NoLegend()
dev.off()

