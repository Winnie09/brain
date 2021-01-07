library(Seurat)
library(here)
setwd(here())
brain <- readRDS('atlasGBM/GBMonly/integrate/harmony/36nonNormalSeuratGene2000/res/humanAtlas_harmony.rds')
clu <- paste0('cluster ', as.character(Idents(brain)))
names(clu) <- colnames(brain@assays$RNA@counts)
clu[clu == 'cluster 5'] <- 'oligodendrocyte'
clu[clu == 'cluster 9'] <- 'microglia'
clu[clu == 'cluster 10'] <- 'endothelial'
clu[clu == 'cluster 13'] <- 'oligodendrocyte'
clu[clu == 'cluster 15'] <- 'oligodendrocyte'
clu.new <- sapply(clu, function(i) ifelse(grepl('cluster', i), NA, i))

atlas.gbm <- readRDS('atlasGBM/humanGBM/integrate/harmony/36nonNormal_allen10xOnly_seuratGene/res/atlasGBM_harmony.rds')
u <- atlas.gbm@reductions$umap@cell.embeddings
meta <- atlas.gbm@meta.data
ct <- meta$celltype
names(ct) <- meta$cell
ct[names(clu.new)] <- clu.new

library(ggplot2)
library(scattermore)
library(RColorBrewer)
pd <- data.frame(umap1 = u[,1], umap2 = u[,2], celltype = ct, feature = sapply(meta$study, function(i) ifelse(grepl('GBM', i), 'GBM', 'atlas')))
  
pd.text <- aggregate(pd[,1:2], list(pd[,3], pd[,4]), median)
colnames(pd.text) <- c('celltype', 'feature','umap1', 'umap2')

p <- ggplot() + geom_scattermore(data  = pd, aes(x=umap1,y=umap2,col=celltype),alpha=0.2,size=0.1) + 
  geom_text(data = pd.text, aes(x = umap1, y = umap2, label = celltype))+
theme_classic() + xlab('UMAP1')+ylab('UMAP2')+
theme(legend.position = 'none',legend.title = element_blank()) + 
guides(color = guide_legend(override.aes = list(size = 4,alpha=1)))+
theme(legend.spacing.x = unit(-0.1, 'cm'),legend.spacing.y = unit(-0.3, 'cm'))+
scale_color_manual(values=colorRampPalette(brewer.pal(9, "Set1"))(length(unique(pd[,3]))+1)[1:length(unique(pd[,3]))]) +
  facet_wrap(~feature)
# ggsave(paste0(pdir,'umap_GBM_atlas.png'),p,height=4, width=5,dpi=200)

