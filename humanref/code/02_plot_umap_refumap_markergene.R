library(Seurat)
ser <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/humanref/res/integrate.rds')
ref <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/humanref/res/reference.rds')

meta <- ser@meta.data

ct <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/integrate/harmony/36nonNormalSeuratGene2000/res/celltype_annotation.rds')  
ser@meta.data <- cbind(meta, celltype=ct[rownames(meta)])

png('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/humanref/plot/alltogether_ref_ser.png', width = 4600, height = 1000, res = 250)
p1 <- DimPlot(ref, reduction = "umap", group.by = "Cell.Type", label = TRUE, label.size = 3, repel = TRUE) 
p2 <- DimPlot(ser, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE, label.size = 3 ,repel = TRUE) 
p3 <- DimPlot(ser, reduction = "ref.umap", group.by = "celltype", label = TRUE, label.size = 3 ,repel = TRUE) 
p1 + p2 + p3
dev.off()


png('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/humanref/plot/GBM_umap_celltype.png', width = 1500, height = 1000, res = 250)
DimPlot(ser, reduction = "umap", group.by = "celltype", label = TRUE, label.size = 3, repel = TRUE) 
dev.off()

png('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/humanref/plot/GBM_refumap_celltype.png', width = 1500, height = 1000, res = 250)
DimPlot(ser, reduction = "ref.umap", group.by = "celltype", label = TRUE, label.size = 3 ,repel = TRUE) 
dev.off()


png('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/humanref/plot/GBM_umap_predictcelltype.png', width = 1500, height = 1000, res = 250)
DimPlot(ser, reduction = "umap", group.by = "predicted.celltype", label = TRUE, label.size = 3, repel = TRUE) 
dev.off()

png('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/humanref/plot/GBM_refumap_predictcelltype.png', width = 1500, height = 1000, res = 250)
DimPlot(ser, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE, label.size = 3 ,repel = TRUE) 
dev.off()


png('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/humanref/plot/reference_newumap_celltype.png', width = 1500, height = 1000, res = 250)
DimPlot(ref, reduction = "umap", group.by = "Cell.Type", label = F, label.size = 3 ,repel = TRUE) #+ NoLegend()
dev.off()



reference <- readRDS('/home-4/zji4@jhu.edu/scratch/GBM/atlas/humanref/data/alltogether_cortex_v3.rds')
meta <- read.table('/home-4/zji4@jhu.edu/scratch/GBM/atlas/humanref/data/allcortex_metadata.tsv',sep='\t',header=T,row.names=1)
for (i in colnames(meta))
  reference@meta.data[,i] <- meta[rownames(reference@meta.data),i]
#Idents(reference) <- 'orig.ident'
png('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/humanref/plot/reference_oldumap_celltype.png', width = 1800, height = 1400, res = 250)
DimPlot(reference, reduction = "umap", group.by = "Cell.Type", label = F, label.size = 3 ,repel = TRUE) 
dev.off()




tb <- read.csv('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/doc/Marker_Gene_Human_Soraia.csv', row.names = 1)
tb <- t(tb)
str(tb)

library(ggplot2)
pdir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/humanref/plot/'
for (i in colnames(tb)){
  print(i)
  png(paste0(pdir, 'GBM_refumap_SB_', i, '.png'), width = round(sqrt(length(tb[,i]))) * 250, height = round(sqrt(length(tb[,i]))) * 250)
  g <- intersect(toupper(tb[,i]), rownames(ser@assays$RNA@counts))
  print(FeaturePlot(ser,  features = g, reduction = 'ref.umap'))
  dev.off()
}
  

for (i in colnames(tb)){
  print(i)
  png(paste0(pdir, 'GBM_umap_SB_', i, '.png'), width = round(sqrt(length(tb[,i]))) * 250, height = round(sqrt(length(tb[,i]))) * 250)
  g <- intersect(toupper(tb[,i]), rownames(ser@assays$RNA@counts))
  print(FeaturePlot(ser,  features = g, reduction = 'umap'))
  dev.off()
}


png(paste0(pdir, 'GBM_NatNeuPaper_marker.png'), width = 1200, height = 1000, res = 200)
FeaturePlot(ser,  features = c('SOX2', 'NES', 'PPP1R17', 'BCL11B'), reduction = 'umap')
dev.off()
png(paste0(pdir, 'GBM_NatNeuPaper_marker_refumap.png'), width = 1200, height = 1000, res = 200)
FeaturePlot(ser,  features = c('SOX2', 'NES', 'PPP1R17', 'BCL11B'), reduction = 'ref.umap')
dev.off()

png(paste0(pdir, 'ref_NatNeuPaper_marker.png'), width = 1200, height = 1000, res = 200)
FeaturePlot(ref,  features = c('SOX2', 'NES', 'PPP1R17', 'BCL11B'), reduction = 'umap')
dev.off()


