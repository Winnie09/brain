key = as.character(commandArgs(trailingOnly = T)[[1]])
# key = 'sub'
# key = 'full'
sp = as.character(commandArgs(trailingOnly = T)[[2]])
setwd(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlas/res/',key))
library(Seurat)
library(ggplot2)
library(cowplot)
u = readRDS('./umap/umap.rds')
meta = readRDS('./meta.rds')
##
clu = readRDS('./cluster/cluster.rds')    
id = setdiff(1:nrow(u),intersect(grep('GBM',meta$study), grep(sp,meta$study,invert=T))) ## find only atlas and sp cells
u = u[id,]
meta = meta[id,]
clu = clu[rownames(u)]

ct = as.character(meta$celltype_super)  

id = intersect(which(!is.na(meta$celltype)), which(ct=='neuronal'))

ct[id] <- as.character(meta$celltype)[id]  ##### detailed neuronal
names(ct) = rownames(u)
ct[grep('GBM',meta$study)] <- sp
ct = factor(ct,levels=c(sp, setdiff(unique(ct),sp)))  ##!!

## color
colv = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/colors3.rds')
colv = colv[1:length(unique(ct))]
names(colv) = c(sp,setdiff(unique(ct),sp))

## clu center
cm = aggregate(u,list(clu),mean)
colnames(cm) <- c('cluster','x','y')
cluct <- sapply(unique(clu),function(i){
  ss = table(ct[names(ct)%in%names(clu[clu==i])])
  ss = ss/sum(ss)
})
colnames(cluct) = unique(clu)
library(reshape2)
pdcluct <- melt(cluct)
colnames(pdcluct) = c('ct','cluster','prop')

##
ctclu <- sapply(unique(ct),function(i){
  ss = table(clu[ct==i])
  v = rep(0,max(clu))
  names(v) <- seq(1,max(clu))
  ss = ss/sum(ss)
  v[names(ss)] = ss
  v
})
colnames(ctclu) = unique(ct)
pdctclu = melt(ctclu)
colnames(pdctclu) = c('cluster','ct','prop')
pdctclu$ct = factor(as.character(pdctclu$ct),levels = levels(pdcluct$ct))
###############################
library(ggplot2)
library(gridExtra)
dir.create(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlas/plot/',key,'/umap_cluster'),showWarnings = F,recursive = T)
p1 <- ggplot() + geom_point(data=data.frame(umap1=u[,1],umap2=u[,2],ct=ct),aes(x=umap1,y=umap2,col=ct),alpha=0.1,size=0.2) + 
  theme_classic() + xlab('UMAP1')+ylab('UMAP2')+
  theme(legend.position = 'right',legend.title = element_blank()) + scale_color_manual(values=colv)+
  guides(color = guide_legend(override.aes = list(size = 3,alpha=1)))+
  theme(legend.spacing.x = unit(-0.1, 'cm'),legend.spacing.y = unit(-0.8, 'cm'), legend.text = element_text(size=7))
p2 <- ggplot() + geom_point(data=data.frame(umap1=u[,1],umap2=u[,2],cluster=as.factor(clu)),aes(x=umap1,y=umap2,col=cluster),alpha=0.5,size=0.2) + 
  geom_text(data=cm,aes(x=x,y=y,label=cluster),size=3) + theme_classic() + xlab('UMAP1')+ylab('UMAP2')+
  theme(legend.position = 'none')
p3 <- ggplot(pdcluct,aes(x=ct,y=cluster,fill=prop)) +
  geom_tile(aes(fill = prop), colour = "black") +
  scale_fill_gradient(low = "white",high = "steelblue")+
  theme_classic() + 
  theme(axis.text.x = element_text(angle=45,hjust=1)) + xlab('') + ylab('cluster') + ggtitle('cluster composition(rowSum=1)')
p4 <- ggplot(pdctclu,aes(x=ct,y=cluster,fill=prop))+
  geom_tile(aes(fill = prop), colour = "black") +
  scale_fill_gradient(low = "white",high = "steelblue")+
  theme_classic() + 
  theme(axis.text.x = element_text(angle=45,hjust=1)) + xlab('') + ylab('cluster') + ggtitle('celltype distribution(colSum=1)')
pdf(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlas/plot/',key,'/umap_cluster/',sp,'.pdf'),width=18,height=12)
grid.arrange(p1,p2,p3,p4,nrow=2)
dev.off()
png(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlas/plot/',key,'/umap_cluster/',sp,'.png'),width=1800,height=1200)
grid.arrange(p1,p2,p3,p4,nrow=2)
dev.off()
