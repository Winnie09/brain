setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/')
source('./resource/function.R')
u = readRDS('./brain/refineSeurat/result/03_cluster_integrate/umap.rds')
u = u[,c(2,1)]
meta = readRDS('./brain/refineSeurat/result/03_cluster_integrate/clu_meta.rds')
fullmeta = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlas/old/20200212/result/full/meta.rds')
meta$celltype = fullmeta[match(meta$cell, fullmeta$cell), 'celltype_super']

plotdir = './brain/refineSeurat/plot/03_cluster_integrate/'


for (i in c("celltype",'study','clu')){
  meta[is.na(meta[,i]),i] = 'unknown'
  tab = table(meta[,i])
  meta[,i] = paste0(meta[,i],'(',tab[match(meta[,i],names(tab))],')')  
}


mat = readRDS('./brain/refineSeurat/result/03_cluster_integrate/mat_c_transform.rds')


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
pdf(paste0(plotdir,'umap_cell_level.pdf'),height=3,width=12)
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
  if (g %in% rownames(mat)){
    print(g)
    plist[[g]] = plotfunc(g)  
  }
}

pdf(paste0(plotdir,'umap_marker_gene_cell_level.pdf'),height=6,width=15)
grid.arrange(grobs=plist,nrow=2)
dev.off()
png(paste0(plotdir,'umap_marker_gene_cell_level.png'),height=600,width=1600)
grid.arrange(grobs=plist,nrow=2)
dev.off()

