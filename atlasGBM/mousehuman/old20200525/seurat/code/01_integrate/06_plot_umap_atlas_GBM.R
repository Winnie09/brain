library(Seurat)
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlas/seurat/')
# key = as.character(commandArgs(trailingOnly = T)[[1]])
key = 'atlas_GBM_homolog_sub'
d = readRDS(paste0('./result/',key,'/integrated.rds')) ## 16482 357321
# fullmeta = readRDS('./old/20200212/result/full/meta.rds')

u = readRDS(paste0('./result/',key,'/umap/umap.rds'))
meta = readRDS(paste0('./result/',key,'/meta.rds'))

plotdir = paste0('./plot/', key,'/')
dir.create(plotdir,showWarnings = F, recursive = T)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(gridExtra)
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
p1 <- ggplot() + geom_point(data=data.frame(umap1=u[,1],umap2=u[,2],study=meta$study_super, atlasGBM = sapply(meta$study_super, function(i) ifelse(grepl('GBM',i),'GBM','Atlas'))),aes(x=umap1,y=umap2,col=study),alpha=0.2,size=0.1) + 
  theme_classic() + xlab('UMAP1')+ylab('UMAP2')+
  theme(legend.position = 'right',legend.title = element_blank()) + 
  guides(color = guide_legend(override.aes = list(size = 4,alpha=1)))+
  theme(legend.spacing.x = unit(-0.1, 'cm'),legend.spacing.y = unit(-0.3, 'cm'))+
   scale_color_manual(values=getPalette(length(unique(meta$study_super))+1)[1:length(unique(meta$study_super))]) + facet_wrap(~atlasGBM)
ggsave(paste0(plotdir,'umap_study.png'),p1,width = 10,height = 4.5,dpi = 200)

id1 = !grepl('unknown',meta$celltype_super) 
p2 <- ggplot() + geom_point(data=data.frame(umap1=u[id1,1],umap2=u[id1,2],study=meta$celltype_super[id1]),aes(x=umap1,y=umap2,col=study),alpha=0.2,size=0.1) + 
  theme_classic() + xlab('UMAP1')+ylab('UMAP2')+
  theme(legend.position = 'right',legend.title = element_blank()) + 
  guides(color = guide_legend(override.aes = list(size = 4,alpha=1)))+
  theme(legend.spacing.x = unit(-0.1, 'cm'),legend.spacing.y = unit(-0.3, 'cm'))+
   scale_color_manual(values=getPalette(length(unique(meta$celltype_super))+1)[1:length(unique(meta$celltype_super))])
ggsave(paste0(plotdir,'umap_celltype_super.png'),p2,width = 9,height = 4.5,dpi = 200)

p3 <- ggplot() + geom_point(data=data.frame(umap1=u[id1,1],umap2=u[id1,2],study=meta$species[id1]),aes(x=umap1,y=umap2,col=study),alpha=0.2,size=0.1) + 
  theme_classic() + xlab('UMAP1')+ylab('UMAP2')+
  theme(legend.position = 'right',legend.title = element_blank()) + 
  guides(color = guide_legend(override.aes = list(size = 4,alpha=1)))+
  theme(legend.spacing.x = unit(-0.1, 'cm'),legend.spacing.y = unit(-0.3, 'cm'))+
   scale_color_manual(values=c('steelblue','orange'))
ggsave(paste0(plotdir,'umap_species.png'),p3,width = 6,height = 4.5,dpi = 200)

id2 = !grepl('unknown',meta$time_super)
p4 <- ggplot() + geom_point(data=data.frame(umap1=u[id2,1],umap2=u[id2,2],study=meta$time_super[id2]),aes(x=umap1,y=umap2,col=study),alpha=0.2,size=0.01) + 
  theme_classic() + xlab('UMAP1')+ylab('UMAP2')+
  theme(legend.position = 'right',legend.title = element_blank()) + 
  guides(color = guide_legend(override.aes = list(size = 4,alpha=1)))+
  theme(legend.spacing.x = unit(-0.1, 'cm'),legend.spacing.y = unit(-0.3, 'cm'))+
   scale_color_manual(values=getPalette(length(unique(meta$time_super))+1)[1:length(unique(meta$time_super))])
ggsave(paste0(plotdir,'umap_time_super.png'),p4, width = 6,height = 4.5,dpi = 200)

pc <- ggplot() + geom_point(data=data.frame(umap1=u[id2,1],umap2=u[id2,2],study=meta$celltype_super[id2]),aes(x=umap1,y=umap2,col=study),alpha=0.2,size=0.1) + 
  theme_classic() + xlab('UMAP1')+ylab('UMAP2')+
  theme(legend.position = 'right',legend.title = element_blank()) + 
  guides(color = guide_legend(override.aes = list(size = 4,alpha=1)))+
  theme(legend.spacing.x = unit(-0.1, 'cm'),legend.spacing.y = unit(-0.3, 'cm'))+
   scale_color_manual(values=getPalette(length(unique(meta$celltype_super))+1)[1:length(unique(meta$celltype_super))])+
  facet_wrap(~study)
ggsave(paste0(plotdir,'umap_celltype.png'),pc, height=6,width=12,dpi=200)


pt <- ggplot() + geom_point(data=data.frame(umap1=u[id2,1],umap2=u[id2,2],study=meta$time_super[id2]),aes(x=umap1,y=umap2,col=study), alpha = 0.1, size = 0.01) + 
  theme_classic() + xlab('UMAP1')+ylab('UMAP2')+
  theme(legend.position = 'right',legend.title = element_blank()) + 
  guides(color = guide_legend(override.aes = list(size = 4,alpha=1)))+
  theme(legend.spacing.x = unit(-0.1, 'cm'),legend.spacing.y = unit(-0.3, 'cm'))+
   scale_color_manual(values=getPalette(length(unique(meta$time_super))+1)[1:length(unique(meta$time_super))])+
  facet_wrap(~study) + scale_alpha_manual(values=c(0.01,0.05,0.05,0.5,0.5))
ggsave(paste0(plotdir,'umap_time.png'), pt, width=7,height=4.5,dpi=200)



ps <- ggplot() + geom_point(data=data.frame(umap1=u[id2,1],umap2=u[id2,2],study=meta$study_super[id2]),aes(x=umap1,y=umap2,col=study),alpha=0.2,size=0.1) + 
  theme_classic() + xlab('UMAP1')+ylab('UMAP2')+
  theme(legend.position = 'right',legend.title = element_blank()) + 
  guides(color = guide_legend(override.aes = list(size = 4,alpha=1)))+
  theme(legend.spacing.x = unit(-0.1, 'cm'),legend.spacing.y = unit(-0.3, 'cm'))+
   scale_color_manual(values=getPalette(length(unique(meta$study_super))+1)[1:length(unique(meta$study_super))])+
  facet_wrap(~study)
ggsave(paste0(plotdir,'umap_study_super.png'), ps, width=10,height=7,dpi=200)

plotfunc <- function(gene){
  p <- ggplot() + geom_point(data=data.frame(umap1=u[,1],umap2=u[,2],expr = mat[gene,]),aes(x=umap1,y=umap2,col=expr),alpha=0.8,size=1) + 
  theme_classic() + xlab('UMAP1')+ylab('UMAP2')+
  ggtitle(gene)+ 
  scale_color_gradientn(colors=brewer.pal(9,'PuBu')) 
  return(p)
}
gvec = c('CD40','CD68','CD11b','CD45',  ## microglia
         'ALDH1L1','GFAP','EAA1','GLAST','EAA2','GLT1','GLT-1', ## astrocyte
         'MBP','OLIG1','OLIG2', ## oligodendrocytes
         'MAP2', 'NeuN', 'PSD95', ## mature neurons
         'SOX10', ## schwann cell precursor, myelinating schwann cell
         'SOX2','HES1', ##  neuroepithelial
         'GAT1', 'GAD65','GAD67', ## GABAergic neurons
         'TH','DAT', ## dopaminnergic neuron
         'TPH',
         'ChAT',
         'POU5F1', 'SOX2', 'KLF4', 'NANOG','MYC','SALL4') ## embryonic
mat = as.matrix(d@assays$RNA@counts) ## count before integration
rm(d)
gvec = intersect(gvec, rownames(mat))
plist = list()
for (g in gvec){
  plist[[g]] = plotfunc(g)
}
ggsave(paste0(plotdir,'umap_marker_gene.png'),grid.arrange(grobs=plist,nrow=2),height=5,width=11,dpi=200)

