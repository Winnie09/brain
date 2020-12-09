#### atlas all umap
umap = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlas/harmony/result/atlasOnly_homolog/umap.rds')
pd = data.frame(umap1 = umap[,1], umap2 =umap[,2], study = sub(':.*','', rownames(umap)))
library(ggplot2)
p <- ggplot() + geom_point(data = pd, aes(x = umap1, y = umap2, color=study), size=0.2, alpha=0.2) + theme_classic() + xlab('UMAP1') + ylab('UMAP2')  + guides(color = guide_legend(override.aes = list(size = 4,alpha=1)))
plotdir = '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlas/harmony/plot/atlasOnly_homolog/'
ggsave(paste0(plotdir,'umap_atlas.png'),
       p,
       height=7,width=9,
       dpi=200)


##### atlas sub umap
res <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlas/harmony/result/atlasOnly_homolog_sub/harmony_embeddings.rds')
library(umap)
set.seed(12345)
umap <- umap(res)$layout
saveRDS(umap, '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlas/harmony/result/atlasOnly_homolog_sub/umap.rds')

meta <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlas/seurat/result/atlasOnly_homolog_sub/meta.rds')
meta = meta[match(rownames(umap), as.character(meta$cell)), ]
pd = data.frame(umap1 = umap[,1], umap2 =umap[,2], study = sub(':.*','', rownames(umap)))
library(ggplot2)
library(gridExtra)
library(cowplot)
library(RColorBrewer)
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
p1 <- ggplot() + geom_point(data = pd, aes(x = umap1, y = umap2, color=study), size=0.2, alpha=0.2) + theme_classic() + xlab('UMAP1') + ylab('UMAP2') + guides(color = guide_legend(override.aes = list(size = 4,alpha=1))) +
  scale_color_manual(values=getPalette(length(unique(meta$study))+1)[1:length(unique(meta$study))])
plotdir = '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlas/harmony/plot/atlasOnly_homolog_sub/'
ggsave(paste0(plotdir,'umap.png'),
       p1,
       height=7,width=9,
       dpi=200)
u = umap
id1 = !grepl('unknown',meta$celltype_super) 
p2 <- ggplot() + geom_point(data=data.frame(umap1=u[id1,1],umap2=u[id1,2],study=meta$celltype_super[id1]),aes(x=umap1,y=umap2,col=study),alpha=0.2,size=0.1) + 
  theme_classic() + xlab('UMAP1')+ylab('UMAP2')+
  theme(legend.position = 'right',legend.title = element_blank()) + 
  guides(color = guide_legend(override.aes = list(size = 4,alpha=1)))+
  theme(legend.spacing.x = unit(-0.1, 'cm'),legend.spacing.y = unit(-0.3, 'cm'))+
   scale_color_manual(values=getPalette(length(unique(meta$celltype_super))+1)[1:length(unique(meta$celltype_super))])
p3 <- ggplot() + geom_point(data=data.frame(umap1=u[id1,1],umap2=u[id1,2],study=meta$species[id1]),aes(x=umap1,y=umap2,col=study),alpha=0.2,size=0.1) + 
  theme_classic() + xlab('UMAP1')+ylab('UMAP2')+
  theme(legend.position = 'right',legend.title = element_blank()) + 
  guides(color = guide_legend(override.aes = list(size = 4,alpha=1)))+
  theme(legend.spacing.x = unit(-0.1, 'cm'),legend.spacing.y = unit(-0.3, 'cm'))+
   scale_color_manual(values=c('orange','steelblue'))

id2 = !grepl('unknown',meta$time_super)
p4 <- ggplot() + geom_point(data=data.frame(umap1=u[id2,1],umap2=u[id2,2],study=meta$time_super[id2]),aes(x=umap1,y=umap2,col=study),alpha=0.2,size=0.1) + 
  theme_classic() + xlab('UMAP1')+ylab('UMAP2')+
  theme(legend.position = 'right',legend.title = element_blank()) + 
  guides(color = guide_legend(override.aes = list(size = 4,alpha=1)))+
  theme(legend.spacing.x = unit(-0.1, 'cm'),legend.spacing.y = unit(-0.3, 'cm'))+
   scale_color_manual(values=getPalette(length(unique(meta$time_super))+1)[1:length(unique(meta$time_super))])
ggsave(paste0(plotdir,'umap_meta.png'),
       grid.arrange(p1,p2,p3,p4,nrow=2, layout_matrix = rbind(c(1,1,2,2,2),c(3,3,4,4,5))),
       height=7,width=17.5,
       dpi=200)


##### atlas all pca
res <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlas/harmony/result/atlasOnly_homolog/pca_embeddings.rds')
varvec <- apply(res, 2, var)
pd = data.frame(pca1 = res[,1], pca2 =res[,2], study = sub(':.*','', rownames(res)))
library(ggplot2)
p <- ggplot() + geom_point(data = pd, aes(x = pca1, y = pca2, color=study), size=0.2, alpha=0.2) + theme_classic() + xlab(paste0('PC1(', round(varvec[1]/sum(varvec),3)*100, '%)')) + ylab(paste0('PC2(', round(varvec[2]/sum(varvec),3)*100, '%)')) + guides(color = guide_legend(override.aes = list(size = 4,alpha=1)))
plotdir = '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlas/harmony/plot/'
ggsave(paste0(plotdir,'pca_atlas.png'),
       p,
       height=7,width=9,
       dpi=200)

##### atlas sub pca
res <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlas/harmony/result/atlasOnly_homolog_sub/pca_embeddings.rds')
varvec <- apply(res, 2, var)
pd = data.frame(pca1 = res[,1], pca2 =res[,2], study = sub(':.*','', rownames(res)))
library(ggplot2)
p <- ggplot() + geom_point(data = pd, aes(x = pca1, y = pca2, color=study), size=0.2, alpha=0.2) + theme_classic() + xlab(paste0('PC1(', round(varvec[1]/sum(varvec),3)*100, '%)')) + ylab(paste0('PC2(', round(varvec[2]/sum(varvec),3)*100, '%)')) + guides(color = guide_legend(override.aes = list(size = 4,alpha=1)))
plotdir = '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlas/harmony/plot/atlasOnly_homolog_sub/'
ggsave(paste0(plotdir,'pca_atlas.png'),
       p,
       height=7,width=9,
       dpi=200)

