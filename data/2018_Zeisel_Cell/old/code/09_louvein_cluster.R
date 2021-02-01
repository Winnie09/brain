setwd("/scratch/users/whou10@jhu.edu/Wenpin/brain/data/2018_Zeisel_Cell/")
library(ggplot2)
umap <- readRDS('./result/umap.rds')
pr = readRDS('./result/pca.rds')
pr = pr$x[,1:20]
suppressMessages(library(scran))
suppressMessages(library(igraph))
graph <- buildSNNGraph(pr,k=20,d=NA,transposed=T)
res <- cluster_louvain(graph)$membership
saveRDS(res,file='./result/cluster.rds')
cm <- aggregate(umap,list(res),mean)
colnames(cm) <- c('cluster','x','y')
pdf('./plot/umap_cluster.pdf',width=6,height=6)
ggplot() + geom_point(data=data.frame(x=umap[,1],y=umap[,2],cluster=as.factor(res)),aes(x=x,y=y,col=cluster),alpha=.2) + geom_text(data=cm,aes(x=x,y=y,label=cluster)) + theme_classic() + xlab('UMAP1')+ylab('UMAP2')+
  theme(legend.position = 'none')
dev.off()
