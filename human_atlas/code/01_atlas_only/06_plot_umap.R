library(Seurat)
ref.integrated <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/human_atlas/result/integrated.rds')
u = ref.integrated@reductions$umap@cell.embeddings
str(u)

meta = ref.integrated@meta.data
u = u[meta$cell,]
## correct
# meta$celltype[meta$celltype == 'oligodendrocytes precursor cells'] <- 'oligodendrocyte precursor cells'
# 
# meta$celltype_super <- sapply(meta$celltype,function(i){
#   if (grepl('embryonic', i)){
#     'embryonic'
#   } else if (grepl('iPSC',i)){
#     'iPSC'
#   } else if (grepl('neu',i) | grepl('neural',i)){
#     'neuronal'
#   } else if (grepl('progenitor',i)){
#     'progenitor'
#   } else if (grepl('unknown',i)){
#     'unknown'
#   } else {
#     i
#   }
# })
# 
# 
# meta$location_super <- sapply(meta$location,function(i){
#   if (grepl('Frontal',i)){
#     'frontal'
#   } else if (grepl('Inferior',i)){
#     'inferior'
#   } else if (grepl('Parietal', i)){
#     'parietal'
#   } else if (grepl('Temporal',i)){
#     'temporal'
#   } else if (grepl('unkown',i)){
#     'unknown'
#   } else {
#     i
#   }
# })
# 
# meta$time_super <- sapply(meta$time,function(i){
#   v  = c('post','embryo','adult','hESC','iPSC','unknown')
#   for (sv in v){
#     if (grepl(sv,i))  return(sv)
#   }
# })
# 
# 
# meta$study = sub(':.*','',meta$cell)
# ref.integrated@meta.data <- meta
# saveRDS(ref.integrated,'/home-4/whou10@jhu.edu/scratch/Wenpin/brain/human_atlas/result/integrated.rds')
###

tmp <- sapply(5:ncol(meta), function(i){
  tab <- table(meta[,i])
  paste0(meta[,i], '(', tab[match(meta[,i],names(tab))],')')
})
meta[,5:ncol(meta)] <- tmp

library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(gridExtra)
getPalette = colorRampPalette(brewer.pal(9, "Set1"))


png('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/human_atlas/plot/atlas/reference_umap.png',height=1000,width=1500)
p1 <- ggplot() + geom_point(data=data.frame(umap1=u[,1],umap2=u[,2],ct=meta$celltype),aes(x=umap1,y=umap2,col=ct),alpha=0.5,size=0.2) + 
  theme_classic() + xlab('UMAP1')+ylab('UMAP2')+
  theme(legend.position = 'right',legend.title = element_blank()) + #scale_color_manual(values=colv)+
  guides(color = guide_legend(override.aes = list(size = 4,alpha=1)))+
  theme(legend.spacing.x = unit(-0.1, 'cm'),legend.spacing.y = unit(-0.3, 'cm'))+
  scale_color_manual(values=getPalette(length(unique(meta$celltype))+1)[1:length(unique(meta$celltype))])


p2 <- ggplot() + geom_point(data=data.frame(umap1=u[,1],umap2=u[,2],location=meta$location),aes(x=umap1,y=umap2,col=location),alpha=0.5,size=0.2) + 
  theme_classic() + xlab('UMAP1')+ylab('UMAP2')+
  theme(legend.position = 'right',legend.title = element_blank()) + #scale_color_manual(values=colv)+
  guides(color = guide_legend(override.aes = list(size = 4,alpha=1)))+
  theme(legend.spacing.x = unit(-0.1, 'cm'),legend.spacing.y = unit(-0.3, 'cm'))+
  scale_color_manual(values=getPalette(length(unique(meta$location))+1)[1:length(unique(meta$location))])


p3 <- ggplot() + geom_point(data=data.frame(umap1=u[,1],umap2=u[,2],study=meta$study),aes(x=umap1,y=umap2,col=study),alpha=0.5,size=0.2) + 
  theme_classic() + xlab('UMAP1')+ylab('UMAP2')+
  theme(legend.position = 'right',legend.title = element_blank()) + #scale_color_manual(values=colv)+
  guides(color = guide_legend(override.aes = list(size = 4,alpha=1)))+
  theme(legend.spacing.x = unit(-0.1, 'cm'),legend.spacing.y = unit(-0.3, 'cm'))+
  scale_color_manual(values=getPalette(length(unique(meta$study))+1)[1:length(unique(meta$study))])

p4 <- ggplot() + geom_point(data=data.frame(umap1=u[,1],umap2=u[,2],time=meta$time),aes(x=umap1,y=umap2,col=time),alpha=0.5,size=0.2) + 
  theme_classic() + xlab('UMAP1')+ylab('UMAP2')+
  theme(legend.position = 'right',legend.title = element_blank()) + #scale_color_manual(values=colv)+
  guides(color = guide_legend(override.aes = list(size = 4,alpha=1)))+
  theme(legend.spacing.x = unit(-0.1, 'cm'),legend.spacing.y = unit(-0.3, 'cm'))+
  scale_color_manual(values=getPalette(length(unique(meta$time))))
grid.arrange(p1,p2,p3,p4,nrow=2)
dev.off()


### facet wrap
png('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/human_atlas/plot/atlas/reference_umap_celltype.png',height=1000,width=1500)
ggplot() + geom_point(data=data.frame(umap1=u[,1],umap2=u[,2],ct=meta$celltype,ct_super=meta$celltype_super),aes(x=umap1,y=umap2,col=ct),alpha=0.5,size=0.2)+
  theme_classic() + xlab('UMAP1')+ylab('UMAP2')+
  theme(legend.position = 'right',legend.title = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 4,alpha=1)))+
  theme(legend.spacing.x = unit(-0.1, 'cm'),legend.spacing.y = unit(-0.3, 'cm'))+
  facet_wrap(~ct_super) +
  scale_color_manual(values=getPalette(length(unique(meta$celltype))))
dev.off()


png('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/human_atlas/plot/atlas/reference_umap_location.png',height=1000,width=1500)
ggplot() + geom_point(data=data.frame(umap1=u[,1],umap2=u[,2],location=meta$location,location_super=meta$location_super),aes(x=umap1,y=umap2,col=location),alpha=0.5,size=0.2)+
  theme_classic() + xlab('UMAP1')+ylab('UMAP2')+
  theme(legend.position = 'right',legend.title = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 4,alpha=1)))+
  theme(legend.spacing.x = unit(-0.1, 'cm'),legend.spacing.y = unit(-0.3, 'cm'))+
  facet_wrap(~location_super) +
  scale_color_manual(values=getPalette(length(unique(meta$location))))
dev.off()


png('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/human_atlas/plot/atlas/reference_umap_time.png',height=1000,width=1500)
ggplot() + geom_point(data=data.frame(umap1=u[,1],umap2=u[,2],time=meta$time,time_super=meta$time_super),aes(x=umap1,y=umap2,col=time),alpha=0.5,size=0.2)+
  theme_classic() + xlab('UMAP1')+ylab('UMAP2')+
  theme(legend.position = 'right',legend.title = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 4,alpha=1)))+
  theme(legend.spacing.x = unit(-0.1, 'cm'),legend.spacing.y = unit(-0.3, 'cm'))+
  facet_wrap(~time_super)+
  scale_color_manual(values=getPalette(length(unique(meta$time))))
dev.off()


png('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/human_atlas/plot/atlas/reference_umap_study.png',height=1000,width=1500)
ggplot() + geom_point(data=data.frame(umap1=u[,1],umap2=u[,2],study=meta$study),aes(x=umap1,y=umap2,col=study),alpha=0.5,size=0.2)+
  theme_classic() + xlab('UMAP1')+ylab('UMAP2')+
  theme(legend.position = 'right',legend.title = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 4,alpha=1)))+
  theme(legend.spacing.x = unit(-0.1, 'cm'),legend.spacing.y = unit(-0.3, 'cm'))+
  facet_wrap(~study)+
  scale_color_manual(values=getPalette(length(unique(meta$study))+1)[1:length(unique(meta$study))])
dev.off()


### more on cell types
df = data.frame(umap1=u[,1],umap2=u[,2],ct=meta$celltype,ct_super=meta$celltype_super)
plist=list()
for (i in unique(df$ct_super)){
  pd = df[df$ct_super==i,]
  plist[[i]] <- ggplot() + geom_point(data=pd,aes(x=umap1,y=umap2,col=ct),alpha=0.5,size=0.2)+
    theme_classic() + xlab('UMAP1')+ylab('UMAP2')+
    theme(legend.position = 'bottom',legend.title = element_blank()) +
    guides(color = guide_legend(override.aes = list(size = 4,alpha=1)))+
    theme(legend.spacing.x = unit(-0.1, 'cm'),legend.spacing.y = unit(-0.3, 'cm'))+
    scale_color_manual(values=getPalette(length(unique(pd$ct))))
}
png('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/human_atlas/plot/atlas/reference_umap_celltype_details.png',height=1200,width=1500)
grid.arrange(grobs=plist)
dev.off()


png('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/human_atlas/plot/atlas/reference_umap_celltype_neuronal.png',height=900,width=800)
pd = df[grep('neuronal',df$ct_super),]
ggplot() + geom_point(data=pd,aes(x=umap1,y=umap2,col=ct),alpha=0.5,size=0.2)+
  theme_classic() + xlab('UMAP1')+ylab('UMAP2')+
  theme(legend.position = 'bottom',legend.title = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 4,alpha=1)))+
  theme(legend.spacing.x = unit(-0.1, 'cm'),legend.spacing.y = unit(-0.3, 'cm'))+
  scale_color_manual(values=getPalette(length(unique(pd$ct))))
dev.off()


png('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/human_atlas/plot/atlas/reference_umap_celltype_neuronal_details.png',height=800,width=1000)
pd = df[grep('neuronal',df$ct_super),]
ggplot() + geom_point(data=pd,aes(x=umap1,y=umap2,col=ct),alpha=0.5,size=0.2)+
  theme_classic() + xlab('UMAP1')+ylab('UMAP2')+
  theme(legend.position = 'bottom',legend.title = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 4,alpha=1)))+
  theme(legend.spacing.x = unit(-0.1, 'cm'),legend.spacing.y = unit(-0.3, 'cm'))+
  facet_wrap(~ct)+
  scale_color_manual(values=getPalette(length(unique(pd$ct))))
dev.off()


png('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/human_atlas/plot/atlas/reference_umap_celltype_nonNeuronal.png',height=1100,width=1000)
pd = df[!grepl('neuronal',df$ct_super),]
ggplot() + geom_point(data=pd,aes(x=umap1,y=umap2,col=ct),alpha=0.5,size=0.2)+
  theme_classic() + xlab('UMAP1')+ylab('UMAP2')+
  theme(legend.position = 'bottom',legend.title = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 4,alpha=1)))+
  theme(legend.spacing.x = unit(-0.1, 'cm'),legend.spacing.y = unit(-0.3, 'cm'))+
  scale_color_manual(values=getPalette(length(unique(pd$ct))))
dev.off()
