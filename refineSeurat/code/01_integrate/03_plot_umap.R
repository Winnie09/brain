library(Seurat)
library(data.table)
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/')
source('./resource/function.R')
u = readRDS('./brain/refineSeurat/result/01_integrate/umap/umap.rds')
meta = readRDS('./brain/refineSeurat/result/01_integrate/meta.rds')

# datafn = paste0('./atlas/result/',key,'/integrated_GBM_reference.rds')
plotdir = './brain/refineSeurat/plot/01_integrate/'
dir.create(plotdir,showWarnings = F, recursive = T)

for (i in c("celltype","location","species","time","gender","celltype_super","location_super","time_super","study")){
  meta[is.na(meta[,i]),i] = 'unknown'
  tab = table(meta[,i])
  meta[,i] = paste0(meta[,i],'(',tab[match(meta[,i],names(tab))],')')  
} 

library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(gridExtra)
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
p1 <- ggplot() + geom_point(data=data.frame(umap1=u[,1],umap2=u[,2],ct=meta$celltype_super),aes(x=umap1,y=umap2,col=ct),alpha=1,size=0.5) + 
  theme_classic() + xlab('UMAP1')+ylab('UMAP2')+
  theme(legend.position = 'right',legend.title = element_blank()) + 
  guides(color = guide_legend(override.aes = list(size = 4,alpha=1)))+
  theme(legend.spacing.x = unit(-0.1, 'cm'),legend.spacing.y = unit(-0.3, 'cm'))+
  scale_color_manual(values=getPalette(length(unique(meta$celltype_super))+1)[1:length(unique(meta$celltype_super))])


p2 <- ggplot() + geom_point(data=data.frame(umap1=u[,1],umap2=u[,2],location=meta$location_super),aes(x=umap1,y=umap2,col=location),alpha=1,size=0.5) + 
  theme_classic() + xlab('UMAP1')+ylab('UMAP2')+
  theme(legend.position = 'right',legend.title = element_blank()) + 
  guides(color = guide_legend(override.aes = list(size = 4,alpha=1)))+
  theme(legend.spacing.x = unit(-0.1, 'cm'),legend.spacing.y = unit(-0.3, 'cm'))+
  scale_color_manual(values=getPalette(length(unique(meta$location_super))+1)[1:length(unique(meta$location_super))])


p3 <- ggplot() + geom_point(data=data.frame(umap1=u[,1],umap2=u[,2],study=meta$study),aes(x=umap1,y=umap2,col=study),alpha=1,size=0.5) + 
  theme_classic() + xlab('UMAP1')+ylab('UMAP2')+
  theme(legend.position = 'right',legend.title = element_blank()) + 
  guides(color = guide_legend(override.aes = list(size = 4,alpha=1)))+
  theme(legend.spacing.x = unit(-0.1, 'cm'),legend.spacing.y = unit(-0.3, 'cm'))+
  scale_color_manual(values=c('lightblue','orange'))

p4 <- ggplot() + geom_point(data=data.frame(umap1=u[,1],umap2=u[,2],time=meta$time),aes(x=umap1,y=umap2,col=time),alpha=1,size=0.5) + 
  theme_classic() + xlab('UMAP1')+ylab('UMAP2')+
  theme(legend.position = 'right',legend.title = element_blank()) + 
  guides(color = guide_legend(override.aes = list(size = 4,alpha=1)))+
  theme(legend.spacing.x = unit(-0.1, 'cm'),legend.spacing.y = unit(-0.3, 'cm'))+
  scale_color_manual(values=getPalette(length(unique(meta$time))))
png(paste0(plotdir,'umap.png'),height=500,width=1000)
grid.arrange(p1,p2,p3,p4,nrow=2)
dev.off()

pdf(paste0(plotdir,'umap.pdf'),height=6,width=11)
grid.arrange(p1,p2,p3,p4,nrow=2)
dev.off()


