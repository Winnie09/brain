key = as.character(commandArgs(trailingOnly = T)[[1]])
# key = 'atlasOnly_homolog_sub'
# key = 'atlas_GBM_homolog_sub'
setwd(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlas/result/',key))
# library(Seurat)
library(ggplot2)
library(cowplot)
u = readRDS('./umap/umap.rds')
meta = readRDS('./meta.rds')
##
clu = readRDS('./cluster/cluster.rds')    
if (exists("sp")){
  id = union(grep('GBM',meta$study, invert=T), grep(sp,meta$study)) ## find only atlas and sp cells
} else {
  id = grep('GBM',meta$study, invert=T) ## find only atlas 
}
u = u[id,]
meta = meta[id,]
clu = clu[rownames(u)]

ct = as.character(meta$celltype_super)  
id = intersect(which(!is.na(meta$celltype)), which(ct=='neuronal'))
ct[id] <- as.character(meta$celltype)[id]  ##### detailed neuronal
names(ct) = rownames(u)
# if (exists("sp")){
#   ct[grep('GBM',meta$study)] <- sp
#   ct = factor(ct,levels=c(sp, setdiff(unique(ct),sp)))  
# } else {
#   ct = factor(ct)  
# }

tab = table(ct)
for (i in names(tab)){
  ct[ct==i] <- paste0(i, '(', tab[i],')')
}

## color
colv = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/colors3.rds')
colv = colv[1:length(unique(clu))]
names(colv) = unique(clu)


## clu center
cm = aggregate(u,list(clu),mean)
colnames(cm) <- c('cluster','x','y')
cluct <- sapply(unique(clu),function(i){
  ss = table(ct[names(ct)%in%names(clu[clu==i])])
  ss = ss/sum(ss)
  tmp <- rep(0,length(unique(ct)))
  names(tmp) <- unique(ct)
  tmp[names(ss)] <- ss
  tmp
})
colnames(cluct) = unique(clu)
library(reshape2)
pdcluct <- melt(cluct)
colnames(pdcluct) = c('ct','cluster','prop')
pdcluct$ct = factor(as.character(pdcluct$ct), levels= sort(unique(as.character(pdcluct$ct))))

## ct center
cm_ct = aggregate(u, list(ct), mean)
colnames(cm_ct) <- c('ct','x','y')

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
dir.create(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlas/plot/',key,'/umap_cluster',max(clu)),showWarnings = F,recursive = T)

p1 <- ggplot() + geom_point(data=data.frame(umap1=u[,1],umap2=u[,2],cluster=as.factor(clu)),aes(x=umap1,y=umap2,col=cluster),alpha=0.3,size=0.2) +scale_color_manual(values=colv)+ 
  geom_text(data=cm,aes(x=x,y=y,label=cluster),size=3) + theme_classic() + xlab('UMAP1')+ylab('UMAP2')+
  theme(legend.position = 'none') 
p2 <- ggplot(pdcluct,aes(x=ct,y=cluster,fill=prop)) +
  geom_tile(aes(fill = prop), colour = "black") +
  scale_fill_gradient(low = "white",high = "steelblue")+
  theme_classic() + 
  theme(axis.text.x = element_text(angle=45,hjust=1)) + xlab('') + ylab('cluster') + ggtitle('cluster composition(rowSum=1)')
p3 <- ggplot(pdctclu,aes(x=ct,y=cluster,fill=prop))+
  geom_tile(aes(fill = prop), colour = "black") +
  scale_fill_gradient(low = "white",high = "steelblue")+
  theme_classic() + 
  theme(axis.text.x = element_text(angle=45,hjust=1)) + xlab('') + ylab('cluster') + ggtitle('celltype distribution(colSum=1)')
if (!exists('sp')) sp = 'umap_cluster'
ggsave(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlas/plot/',key,'/umap_cluster', max(clu),'/',sp,'.png'),
       grid.arrange(p1,p2,p3,nrow=1),
       width=15,height=4,
       dpi=200)


