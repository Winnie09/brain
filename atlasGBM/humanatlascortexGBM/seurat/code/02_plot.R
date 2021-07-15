library(ggplot2)
d <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/humanatlascortexGBM/seurat/res/project.rds')
p1 <- DimPlot(r, reduction = "umap", group.by = "seurat_clusters", label = TRUE, label.size = 3,
              repel = TRUE) + NoLegend() + ggtitle("Reference annotations")
p2 <- DimPlot(d, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE,
              label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Query transferred labels") + xlim(c(min(r$umap@cell.embeddings[,1]),max(r$umap@cell.embeddings[,1]))) + ylim(c(min(r$umap@cell.embeddings[,2]),max(r$umap@cell.embeddings[,2])))
pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/humanatlascortexGBM/seurat/plot/umap.pdf',width=14,height=7)
print(p1 + p2)
dev.off()
