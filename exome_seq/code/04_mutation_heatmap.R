library(pheatmap)
library(reshape2)
library(ggplot2)
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/')
mut <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/exome_seq/data/mutation.rds')
m <- readRDS('atlasGBM/GBMonly/integrate/seurat/36nonNormalSeuratGene2000/res/meta.rds')
d <- readRDS('atlasGBM/GBMonly/integrate/seurat/36nonNormalSeuratGene2000/res/umap_embeddings.rds')
cl <- readRDS('atlasGBM/GBMonly/integrate/seurat/36nonNormalSeuratGene2000/res/cluster.rds')
ct <- read.csv('doc/GBM_cluster_annotations.csv',as.is=T)

m <- unique(m[,-c(1:4)])
rownames(m) <- m$Sample.ID
m <- m[,c('Location','Pathology','Tumor.Grade','Treatment')]

pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/exome_seq/plot/mutation_heatmap.pdf')
heatmap(mut,scale='none')
dev.off()

mut <- mut[rownames(mut)!='GBM048',]
mut <- mut[,colSums(mut) > 0]
pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/exome_seq/plot/mutation_heatmap_rm048.pdf')
heatmap(mut,scale='none')
dev.off()

