u <- readRDS('Downloads/umap.rds')
m <- readRDS('Downloads/meta.rds')
library(ggplot2)
id <- which(!is.na(m[,'celltype_super']))
ggplot(data.frame(x=u[id,1],y=u[id,2],ct=m[id,'celltype_super']),aes(x=x,y=y,col=ct)) + geom_point() + theme_classic() + facet_wrap(~ct)
library(TSCAN)
clu <- kmeans(u,12,iter.max = 1000)$cluster
em <- exprmclust(t(u),reduce=F,cluster = clu)
nn <- unique(m[,'celltype_super'])
cluct <- sapply(1:max(clu),function(i) {
  tab <- table(m[clu==i,'celltype_super'])
  v <- rep(0,length(nn))
  names(v) <- nn
  v[names(tab)] <- as.vector(tab)
  v
})
row.names(cluct) <- nn
row.names(cluct)[apply(cluct/rowSums(cluct),2,which.max)]
row.names(cluct)[apply(cluct,2,which.max)]

clucenter <- aggregate(u,list(clu),mean)
colnames(clucenter) <- c('cluster','x','y')

library(igraph)
edgelist <- matrix(as.numeric(as_edgelist(em$MSTtree)),ncol=2)
ap <- data.frame(x=clucenter[edgelist[,1],'x'],y=clucenter[edgelist[,1],'y'],xend=clucenter[edgelist[,2],'x'],yend=clucenter[edgelist[,2],'y'])

ggplot() + geom_point(data=data.frame(x=u[id,1],y=u[id,2],ct=m[id,'celltype_super']),aes(x=x,y=y,col=ct)) + geom_text(data=clucenter,aes(x=x,y=y,label=cluster)) + geom_segment(data=ap,aes(x=x,y=y,xend=xend,yend=yend)) + theme_classic()

save.image('Downloads/tscan.rda')

