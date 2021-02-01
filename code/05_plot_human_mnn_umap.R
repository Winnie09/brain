u = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/result/merge_human/mnn_umap.rds')
batch <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/result/merge_human/batch.rds')
batch = rep(names(batch), batch)
library(ggplot2)
pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/plot/merge_human/mnn_umap.pdf',width=6,height=3)
ggplot() + geom_point(data=data.frame(umap1=u[,1],umap2=u[,2],batch=as.factor(batch)),aes(x=umap1,y=umap2,col=batch),alpha=0.2,size=0.1) + theme_classic() + guides(color = guide_legend(override.aes = list(size = 5,alpha=1)))+
  xlim(-9,10)+ylim(-8,8)
dev.off()


