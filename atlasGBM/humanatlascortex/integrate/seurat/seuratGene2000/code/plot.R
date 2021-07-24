library(Seurat)
library(ggplot2)
library(pheatmap)
d <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/humanatlascortex/integrate/seurat/seuratGene2000/res/humanAtlas_harmony.rds')

gl <- read.csv('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/humanatlascortex/integrate/plot/gl/Marker_Gene_Tim.csv',as.is=T)
#gl2 <- read.csv('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/humanatlascortex/plot/gl/Marker_Gene_Human_Soraia.csv',as.is=T)
gl2 <- c('DLX1','DLX2','CALR','NELL1','SOX2','PAX6','EOMES','NEUROD1')


## within each dataset (study) and cluster, scale the data
clu <- Idents(d)
n <- names(clu)
clu <- as.character(clu)
names(clu) <- n
m <- d$RNA@counts
stu <- sub(';.*','',colnames(m))
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

clu <- sub('.*_','',colnames(acm))
accm <- sapply(unique(clu),function(i) {
  rowMeans(acm[,clu==i],na.rm=T)
})
accm <- accm[,order(as.numeric(colnames(accm)))]
accm <- accm[rownames(m)[rowMeans(m > 0) >= 0.01],]
fcgl <- sapply(colnames(accm),function(i) {
  names(head(sort(accm[,i]-rowMeans(accm[,colnames(accm)!=i]),decreasing = T),5))
})

pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/humanatlascortex/integrate/plot/plot/fc.pdf',height=20)
pheatmap(accm[as.vector(fcgl),],cluster_rows = F,cluster_cols = F,gaps_row = 5*(1:ncol(accm)))
dev.off()

## marker set 1
gv <- intersect(rownames(acm),unique(c(gl2,unlist(gl))))
cm <- acm[gv,]
clu <- sub('.*_','',colnames(cm))
ccm <- sapply(unique(clu),function(i) {
  rowMeans(cm[,clu==i],na.rm=T)
})

pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/humanatlascortex/integrate/plot/plot/marker.pdf')
pheatmap(ccm)
dev.off()

## marker set 2
gl3 <- toupper(c('Sfrp1', 'Rbfox1', 'Hopx', 'Fam107a', 'Tnc', 'Nhlh1', 'Nhlh2', 'Tp53i11', 'Sla1c3', 'Ppp1r17', 'Ttf1', 'Dlx1', 'Ptprc', 'P2ry12', 'Olig1', 'Col20a1', 'Pmp2', 'Neurod2', 'Slco1c1', 'Pde4dip', 'Aqp4'))
# gv <- intersect(rownames(acm),unique(c(gl3,unlist(gl))))
gv <- intersect(rownames(acm),gl3)
cm <- acm[gv,]
clu <- sub('.*_','',colnames(cm))
ccm <- sapply(unique(clu),function(i) {
  rowMeans(cm[,clu==i],na.rm=T)
})
pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/humanatlascortex/integrate/seurat/seuratGene2000/plot/marker_set2.pdf')
pheatmap(ccm)
dev.off()


## plot UMAP
pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/humanatlascortex/integrate/plot/plot/cluster_umap.pdf')
DimPlot(d,label=T) + theme(legend.position='none')
dev.off()

stu <- sub(';.*','',colnames(d))
tab <- table(stu,Idents(d))
tab <- tab/rowSums(tab)
pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/humanatlascortex/integrate/plot/plot/cluster_prop.pdf')
pheatmap(tab)
dev.off()


stu <- sub(';.*','',colnames(d))
meta <- sapply(unique(stu),function(i) {
  if (i=='2020_Fan_SciAdvA35') i <- '2020_Fan_SciAdv'
  tmp <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/',i,'/proc/meta.rds'))
  colnames(tmp) <- tolower(colnames(tmp))
  colnames(tmp) <- sub('time','age',colnames(tmp))
  data.frame(cell=as.character(tmp[,'cell']),age=as.character(tmp[,'age']),stringsAsFactors = F)
},simplify = F)
meta <- do.call(rbind,meta)
meta$age[meta$age!='adult'] <- 'dev'

pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/humanatlascortex/integrate/plot/plot/dev_umap.pdf')
ggplot(data.frame(umap1=d$umap@cell.embeddings[,1],umap2=d$umap@cell.embeddings[,2],age=meta[match(sub('.*;','',row.names(d$umap@cell.embeddings)),meta[,1]),'age'],stringsAsFactors = F),aes(x=umap1,y=umap2,col=age)) + geom_point(size=0.1,alpha=0.3) + theme_classic()
dev.off()

tab <- table(Idents(d),meta[match(sub('.*;','',colnames(d)),meta[,1]),'age'])
tab <- tab/rowSums(tab)
pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/humanatlascortex/integrate/plot/plot/dev_prop.pdf')
pheatmap(tab)
dev.off()






