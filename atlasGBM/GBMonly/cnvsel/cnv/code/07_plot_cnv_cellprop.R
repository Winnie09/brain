library(pheatmap)
library(reshape2)
library(ggplot2)
pdir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnvsel/cnv/plot/'
m <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/integrate/seurat/36nonNormalSeuratGene2000/res/meta.rds')
d <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/integrate/seurat/36nonNormalSeuratGene2000/res/umap_embeddings.rds')
cl <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/integrate/seurat/36nonNormalSeuratGene2000/res/cluster.rds')
ct <- read.csv('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/doc/GBM_cluster_annotations.csv',as.is=T)

m <- unique(m[,-c(1:4)])
rownames(m) <- m$Sample.ID
m <- m[,c('Location','Pathology','Tumor.Grade','Treatment')]

identical(names(cl),rownames(d))
s <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnvsel/cnv/summary/summary.rds')

type <- rep('NA',nrow(d))
names(type) <- rownames(d)
type[names(s)] <- s
cl_ct <- ct[cl,3] #
names(cl_ct) <- names(cl)

##### Within each GBM sample: cell type * cnv pattern heatmap
ap = sub('_.*', '', names(cl))
for (p in unique(ap)){
  print(p)
  tmp = table(cl_ct[ap == p], type[ap == p])
  tmp <- tmp/sum(ap == p)
  tmp <- tmp[, colSums(tmp) > 0.05, drop = FALSE]
  if (ncol(tmp) >=1 & nrow(tmp) >= 1){
    tmp <- tmp[sort(rownames(tmp)), sort(colnames(tmp))]
    pdf(paste0(pdir, p, '_cnv_ct_cellprop.pdf'), width = 5, height = 5)
    pheatmap(tmp, cluster_cols = F, cluster_rows = F)
    dev.off()
  }
}

##### compare exome-seq to infered cnv profiles: sample * cnv heatmap
exome_mat <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/exome_seq/plot/exome_greater_than_50_percent.rds')

cnv_mat <- matrix(0, nrow = nrow(exome_mat), ncol = ncol(exome_mat))
dimnames(cnv_mat) <- dimnames(exome_mat)

ap = sub('_.*', '', names(cl))
for (p in rownames(cnv_mat)){
  tmpv = sapply(colnames(cnv_mat), function(i)  sum(grepl(i, type[ap == p])))
  cnv_mat[p, ] <- tmpv/sum(tmpv)
}
  
pdf(paste0(pdir, 'infered_cnv_cell_prop.pdf'), height = 4, width = 7.2)
pheatmap(cnv_mat, cluster_cols = F, cluster_rows = F)
dev.off()

