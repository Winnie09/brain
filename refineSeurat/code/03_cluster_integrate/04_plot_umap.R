library(Seurat)
library(data.table)
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/')
source('./resource/function.R')
u = readRDS('./brain/refineSeurat/result/03_cluster_integrate/umap/umap.rds')
meta = readRDS('./brain/refineSeurat/result/03_cluster_integrate/meta.rds')
meta$study = sub(':.*','',meta$cell)
plotdir = './brain/refineSeurat/plot/03_cluster_integrate/'
dir.create(plotdir,showWarnings = F, recursive = T)

for (i in c("celltype")){
  meta[is.na(meta[,i]),i] = 'unknown'
  tab = table(meta[,i])
  meta[,i] = paste0(meta[,i],'(',tab[match(meta[,i],names(tab))],')')  
} 
dlist <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/refineSeurat/result/03_cluster_integrate/dlist.rds')
mat1 = as.matrix(dlist[['2018_Fan_CellResearch']]@assays$RNA@counts)
mat2 = as.matrix(dlist[['2016_La']]@assays$RNA@counts)
mat = cbind(mat1,mat2)


library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(gridExtra)
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
p1 <- ggplot() + geom_point(data=data.frame(umap1=u[,1],umap2=u[,2],ct=meta$celltype),aes(x=umap1,y=umap2,col=ct),alpha=1,size=0.5) + 
  theme_classic() + xlab('UMAP1')+ylab('UMAP2')+
  theme(legend.position = 'right',legend.title = element_blank()) + 
  guides(color = guide_legend(override.aes = list(size = 4,alpha=1)))+
  theme(legend.spacing.x = unit(-0.1, 'cm'),legend.spacing.y = unit(-0.3, 'cm'))+
  scale_color_manual(values=getPalette(length(unique(meta$celltype))+1)[1:length(unique(meta$celltype))])


p2 <- ggplot() + geom_point(data=data.frame(umap1=u[,1],umap2=u[,2],study=meta$study),aes(x=umap1,y=umap2,col=study),alpha=1,size=0.5) + 
  theme_classic() + xlab('UMAP1')+ylab('UMAP2')+
  theme(legend.position = 'right',legend.title = element_blank()) + 
  guides(color = guide_legend(override.aes = list(size = 4,alpha=1)))+
  theme(legend.spacing.x = unit(-0.1, 'cm'),legend.spacing.y = unit(-0.3, 'cm'))+
  scale_color_manual(values=c('lightblue','orange'))
pdf(paste0(plotdir,'umap.pdf'),height=3,width=9)
grid.arrange(p1,p2,layout_matrix=rbind(c(1,1,1,2,2,2), c(1,1,1,2,2,2)))
dev.off()

plotfunc <- function(gene){
  p <- ggplot() + geom_point(data=data.frame(umap1=u[,1],umap2=u[,2],expr = mat[gene,]),aes(x=umap1,y=umap2,col=expr),alpha=1,size=0.5) + 
  theme_classic() + xlab('UMAP1')+ylab('UMAP2')+
  ggtitle(gene)+ 
  scale_color_gradientn(colors=brewer.pal(9,'Blues')) 
  return(p)
}

plist = list()
for (g in c('CD40','CD68','ALDH1L1','MBP','OLIG1','OLIG2','OLIG3','MAP2','SOX10','SOX2')){
  plist[[g]] = plotfunc(g)
}

pdf(paste0(plotdir,'umap_marker_gene.pdf'),height=6,width=17)
grid.arrange(grobs=plist,nrow=2)
dev.off()
png(paste0(plotdir,'umap_marker_gene.png'),height=600,width=1800)
grid.arrange(grobs=plist,nrow=2)
dev.off()

