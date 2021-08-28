library(ggplot2)
library(Seurat)
d = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/humanatlascortexGBM/seurat/res/project_with_celltype_anno.rds')
r = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/humanatlascortex/integrate/seurat/seuratGene2000/res/humanAtlas_harmony_celltype_anno.rds')
p1 <- DimPlot(r, reduction = "umap", group.by = "celltype.anno", label = TRUE, label.size = 3,
              repel = TRUE) + NoLegend() + ggtitle("Human cortex atlas")


p2 <- DimPlot(d, reduction = "ref.umap", group.by = "predicted.celltype_anno", label = TRUE,
              label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Projected GBM cells") + xlim(c(min(r$umap@cell.embeddings[,1]),max(r$umap@cell.embeddings[,1]))) + ylim(c(min(r$umap@cell.embeddings[,2]),max(r$umap@cell.embeddings[,2])))
pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/humanatlascortexGBM/seurat/plot/umap_celltype.pdf',width=12,height=5)
print(p1 + p2)
dev.off()


