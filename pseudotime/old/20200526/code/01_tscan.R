setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlas/result/atlas_GBM_homolog_sub/')
# setwd('/Users/wenpinhou/Dropbox/brain/atlas/result/sub/')
# setwd('/Users/wenpinhou/Dropbox/brain/atlas/old/20200212/result/sub/')
u <- readRDS('./umap/umap.rds')
m <- readRDS('./meta.rds')
m = m[match(rownames(u), m$cell),]
library(ggplot2)
id <- which(!is.na(m[,'celltype_super']))
clu = readRDS('./cluster/cluster.rds')
# ggplot(data.frame(x=u[id,1],y=u[id,2],ct=m[id,'celltype_super']),aes(x=x,y=y,col=ct)) + geom_point() + theme_classic() + facet_wrap(~ct)
library(TSCAN)
set.seed(12345)
# clu <- kmeans(u,12,iter.max = 1000)$cluster

em <- exprmclust(t(u),reduce=F,cluster = clu)
nn <- unique(m[,'celltype_super'])
cluct <- sapply(1:max(clu),function(i) {
  tmp = m[clu==i,'celltype_super']
  tmp = tmp[!grepl('unknown',tmp)]
  tab <- table(tmp)
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
cluct['embryonic',]

library(igraph)
edgelist <- matrix(as.numeric(as_edgelist(em$MSTtree)),ncol=2)
ap <- data.frame(x=clucenter[edgelist[,1],'x'],y=clucenter[edgelist[,1],'y'],xend=clucenter[edgelist[,2],'x'],yend=clucenter[edgelist[,2],'y'])

pd1 = data.frame(x=u[id,1],y=u[id,2],ct=m[id,'celltype_super'])
saveRDS(list(pd1,ap),'/scratch/users/whou10@jhu.edu/Wenpin/brain/pseudotime/plot/plot/pseudotime_umap_pd.rds')
# pdf('/scratch/users/whou10@jhu.edu/Wenpin/brain/pseudotime/plot/plot/pseudotime.pdf',width=10,height=10)
p<- ggplot() + geom_point(data=pd1,aes(x=x,y=y,col=ct),alpha=0.5,size=0.2) + 
  geom_segment(data=ap,aes(x=x,y=y,xend=xend,yend=yend),size=1,color='brown') + 
  geom_text(data=clucenter,aes(x=x,y=y,label=cluster),size=7) + 
  guides(color = guide_legend(override.aes = list(size = 4,alpha=1)))+
  theme_classic() + theme(legend.title = element_blank()) 
# dev.off()
ggsave('/scratch/users/whou10@jhu.edu/Wenpin/brain/pseudotime/plot/plot/pseudotime.png',
       p,
       width=10,height=10,
       dpi=200)

# 
# # load('Downloads/tscan.rda')
# order <- list()
# order[['6_4']] <- TSCANorder(em,MSTorder = c(6,10,4),orderonly = T)
# order[['6_3']] <- TSCANorder(em,MSTorder = c(6,10,3),orderonly = T)
# order[['6_11']] <- TSCANorder(em,MSTorder = c(6,8,11),orderonly = T)
# order[['6_12']] <- TSCANorder(em,MSTorder = c(6,1,12),orderonly = T)
# order[['6_2']] <- TSCANorder(em,MSTorder = c(6,1,2),orderonly = T)
# order[['6_9']] <- TSCANorder(em,MSTorder = c(6,1,7,9),orderonly = T)
# order[['6_5']] <- TSCANorder(em,MSTorder = c(6,1,7,5),orderonly = T)
# 
# 
# meta2 = read.csv('/home-4/whou10@jhu.edu/scratch/Wenpin/GBM/doc/Meta_Wenpin_10292019.csv') 
# meta2$Sample.ID = sapply(meta2$Sample.ID, function(i) sub(' ','',i))
# meta2 = meta2[,!colnames(meta2)%in%c('Progression.date')]
# meta2 = meta2[order(meta2$Location),]
# up <- unique(sub('_.*','',grep('GBM',row.names(u),value=T)))
# pos <- lapply(order,function(s) {
#   v <- 1:length(s)
#   sp <- sub('_.*','',s)
#   tmp <- NULL
#   for (ssp in up) {
#     tmp <- rbind(tmp,data.frame(patient=ssp,pseudotime=v[sp==ssp],stringsAsFactors = F))
#   }
#   df = cbind(tmp, meta2[match(tmp$patient, meta2$Sample.ID),])
#   df$patient = factor(df$patient,levels=meta2[meta2$Pathology=='GBM','Sample.ID'])
#   df
# })
# names(pos) <- names(order)
# saveRDS(pos, './pseudotime/plot/plot/pseudotime_order.rds')
# 
# pdf('./pseudotime/plot/plot/ptdist.pdf',width=4,height=9)
# for (i in names(pos))
# print(ggplot(pos[[i]],aes(x=pseudotime,fill=Location),alpha=0.5) + geom_histogram() + theme_classic() + facet_wrap(~patient,scales = 'free_y',ncol=1,strip.position="right") + theme(strip.text.y = element_text(angle = 0)) + ggtitle(i)+
#         theme(axis.text=element_text(size=5))+
#         scale_fill_brewer(palette='Set2'))
# dev.off()
# 

