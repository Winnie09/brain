library(Seurat)
library(ggplot2)
library(pheatmap)
d <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/integrate/seurat/36nonNormalSeuratGene2000/res/humanAtlas_harmony.rds')

n = names(d@active.ident)
d@active.ident <- factor(as.numeric(d@active.ident), levels = seq(1, length(unique(d@active.ident))))
names(d@active.ident) = n

gl <- read.csv('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/humanatlascortex/integrate/gl/Marker_Gene_Tim.csv',as.is=T)
#gl2 <- read.csv('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/humanatlascortex/plot/gl/Marker_Gene_Human_Soraia.csv',as.is=T)

gl2 <- c('DLX1','DLX2','CALR','NELL1','SOX2','PAX6','EOMES','NEUROD1')

clu <- Idents(d)
table(clu)

n <- names(clu)
clu <- as.character(clu)
names(clu) <- n
m <- d$RNA@counts
stu <- sub('_.*','',colnames(m))

library(Matrix)
acm <- sapply(unique(stu),function(i) {
  tmp <- m[,stu==i]
  tmpclu <- clu[colnames(tmp)]
  cm <- sapply(unique(tmpclu),function(j) {
    rowMeans(tmp[,tmpclu==j,drop=F])
  })
  cm <- t(apply(cm,1,scale))
  colnames(cm) <- paste0(i,'_',unique(tmpclu))
  cm
})
acm <- do.call(cbind,acm)
saveRDS(acm, '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/integrate/seurat/36nonNormalSeuratGene2000/res/acm.rds')

clu <- sub('.*_','',colnames(acm))
accm <- sapply(unique(clu),function(i) {
  rowMeans(acm[,clu==i],na.rm=T)
})
accm <- accm[,order(as.numeric(colnames(accm)))]
accm <- accm[rownames(m)[rowMeans(m > 0) >= 0.01],]
fcgl <- sapply(colnames(accm),function(i) {
  names(head(sort(accm[,i]-rowMeans(accm[,colnames(accm)!=i]),decreasing = T),5))
})

cluanno = read.csv('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/doc/GBM_cluster_annotations.csv', as.is = T)
colnames(accm) = cluanno[,3]

pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/integrate/seurat/36nonNormalSeuratGene2000/plot/fc.pdf',height=20)
pheatmap(accm[as.vector(fcgl),],cluster_rows = F,cluster_cols = F,gaps_row = 5*(1:ncol(accm)))
dev.off()


pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/integrate/seurat/36nonNormalSeuratGene2000/plot/fc.pdf',width=23, height = 6)
pheatmap(t(accm[as.vector(fcgl),]),cluster_rows = F,cluster_cols = F,gaps_col = 5*(1:ncol(accm)))
dev.off()

#######  Here is the new gene list from Tim, need to have their expression heatmap
genelist = readLines('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/doc/GBM_Tim_genes.txt')
genelist = toupper(genelist)
genelist = genelist[genelist %in% rownames(accm)]
pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/integrate/seurat/36nonNormalSeuratGene2000/plot/Hm_GBM_genes_from_Tim.pdf',height=20)
pdf('Hm_GBM_genes_from_Tim.pdf',height=12)
pheatmap(accm[as.vector(genelist),],cluster_rows = F,cluster_cols = F)
dev.off()

#############################
gv <- intersect(rownames(acm),unique(c(gl2,unlist(gl))))
cm <- acm[gv,]
clu <- sub('.*_','',colnames(cm))
ccm <- sapply(unique(clu),function(i) {
  rowMeans(cm[,clu==i],na.rm=T)
})

pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/integrate/seurat/36nonNormalSeuratGene2000/plot/marker.pdf')
pheatmap(ccm)
dev.off()

pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/integrate/seurat/36nonNormalSeuratGene2000/plot/cluster_umap.pdf')
DimPlot(d,label=T) + theme(legend.position='none')
dev.off()

stu <- sub('_.*','',colnames(d))
tab <- table(stu,Idents(d))
tab <- tab/rowSums(tab)
pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/integrate/seurat/36nonNormalSeuratGene2000/plot/cluster_prop.pdf')
pheatmap(tab)
dev.off()



