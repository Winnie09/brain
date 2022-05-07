clu <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/integrate/seurat/36nonNormalSeuratGene2000/res/cluster.rds')
anno <- read.csv('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/doc/GBM_cluster_annotations.csv')
ct <- anno[match(clu, anno[,1]),3]
names(ct) <- names(clu)
af <- list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_myelymreference/cnvclu')
pdir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_myelymreference/compare_cnvclu_with_GBMclu/'
library(ggplot2)
for (f in af) {
  cnv <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_myelymreference/cnvclu/',f))
  df <- data.frame(table(cnv,ct[names(cnv)]))
  colnames(df) <- c('cnv','expr','count')
  df[,1] <- as.character(df[,1])
  tab = table(cnv)
  df[,1] <- paste0(df[,1],'(', tab[df[,1]],')')
  df[,1] <- factor(df[,1], levels = unique(df[,1]))
  
  pdf(paste0(pdir, sub('.rds','.pdf', f)), width = ifelse(f %in% c('GBM074.rds'), 5,4), height=4)
  print(ggplot(df,aes(x=cnv,y=expr,fill=count)) + geom_tile() + theme_classic() + scale_fill_gradient(low='black',high='yellow') + xlab('cnv clusters') + ylab('expression clusters') + theme(axis.text.x = element_text(angle=45, hjust = 1)))
  dev.off()
}

