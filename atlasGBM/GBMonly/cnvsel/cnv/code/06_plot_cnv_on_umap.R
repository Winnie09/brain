library(ggplot2)
library(cowplot)
library(scattermore)
library(RColorBrewer)
library(ggrepel) 

source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/ggplot_theme.R')
theme_set(.new_theme)

u = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/integrate/seurat/36nonNormalSeuratGene2000/res/umap_embeddings.rds')
s = readRDS('/home-4/zji4@jhu.edu/scratch/cnvsel/cnv/summary/summary.rds')


## umap of all samples, with/without cnv, colored by orange/blue
cnv_bi = sapply(as.character(rownames(u) %in% names(s)), function(i) ifelse(i == 'TRUE', 'CNV', 'no_CNV'))
pd <- data.frame(umap1 = u[,1], umap2 = u[,2], cnv = cnv_bi, stringsAsFactors = F)
pd.text <- aggregate(pd[,1:2], list(pd[,3]), median)
colnames(pd.text) <- c('cnv','umap1', 'umap2')
colv = c('royalblue', 'orange')
names(colv) = c('no_CNV', 'CNV')

p <- ggplot() + geom_scattermore(data  = pd, aes(x=umap1,y=umap2,col=cnv),alpha=0.2,size=0.1) + 
  geom_text_repel(data = pd.text, aes(x = umap1, y = umap2, label = cnv))+
theme_classic() + xlab('UMAP1')+ylab('UMAP2')+
theme(legend.position = 'right',legend.title = element_blank()) + 
guides(color = guide_legend(override.aes = list(size = 4,alpha=1)))+
theme(legend.spacing.x = unit(-0.1, 'cm'),legend.spacing.y = unit(-0.3, 'cm'))+
scale_color_manual(values=colv[unique(pd[,3])]) 
ggsave('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnvsel/cnv/plot/cnv_on_umap_TF.pdf', p, height=3, width=3.5,dpi=200)

p <- p + facet_wrap(~cnv)
ggsave('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnvsel/cnv/plot/cnv_on_umap_TF_facet.pdf', p, height=2.6, width=5.5,dpi=200)

## umap of each cnv pattern, use all samples' cells, colored by red/grey
length(unique(s))
# [1] 17872
head(sort(table(s), decreasing = T), 20)
       #       chr10-        chr10-;chr7+       chr10-;chr13- chr10-;chr13-;chr7+ 
       #        20050               14898                9003                5063 
       #       chr21+              chr22-        chr10-;chr6-              chr13- 
       #         4786                4359                3229                3171 
       # chr1-;chr19-               chr7+ chr10-;chr22-;chr7+               chr1- 
       #         2786                2595                2332                2223 
       #        chr6-              chr18+       chr10-;chr22-              chr19- 
       #         2127                1893                1768                1555 
       # chr22-;chr7+              chr11+  chr10-;chr6-;chr7+ chr10-;chr20+;chr7+ 
       #         1471                1372                1352                1336 
sum(sort(table(s), decreasing = T)> 1000)
# [1] 26

large_clone <- names(sort(table(s)[table(s) > 1000], decreasing = TRUE))
s_large_clone <- s[s %in% large_clone]
tb = table(s_large_clone)
pd <- data.frame(cnv = names(tb), numCell = as.vector(tb), stringsAsFactors = F)
pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnvsel/cnv/plot/cnv_larger_than_1k_cells_hist.pdf', width = 3.8, height = 3)
ggplot(data = pd, aes(x = cnv, y = numCell)) + geom_bar(stat = 'identity') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ggtitle(paste0('26 CNV w/ >= 1k cells, in total ', length(s_large_clone), ' cells (out of ', length(s), ' w/ CNV)'))
dev.off()

## plot each CNV pattern on UMAP
pd <- data.frame(umap1 = u[names(s),1], umap2 = u[names(s),2], cnv = s, stringsAsFactors = F)
pd <- pd[pd[,'cnv'] %in%  large_clone, ]

colv <- colorRampPalette(brewer.pal(9, "Set1"))(length(unique(pd[,3]))+1)[1:length(unique(pd[,3]))]
names(colv) <- unique(pd[,3])
p <- ggplot() + geom_scattermore(data  = pd, aes(x=umap1,y=umap2,col=cnv),alpha=1,size=0.8) + 
theme_classic() + xlab('UMAP1')+ylab('UMAP2')+
theme(legend.position = 'right',legend.title = element_blank()) + 
guides(color = guide_legend(override.aes = list(size = 4,alpha=1)))+
theme(legend.spacing.x = unit(-0.1, 'cm'),legend.spacing.y = unit(-0.3, 'cm'))+
scale_color_manual(values=colv[unique(pd[,3])]) 
ggsave('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnvsel/cnv/plot/cnv_on_umap.pdf', p, height=4, width=7,dpi=300)

p <- p + facet_wrap(~cnv)
ggsave('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnvsel/cnv/plot/cnv_on_umap_facet.pdf', p, height=9, width=12,dpi=200)

