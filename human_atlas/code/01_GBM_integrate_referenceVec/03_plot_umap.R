library(Seurat)
datafn = '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/human_atlas/result/referenceVec/GBM_integratedUMAP.rds'
plotdir = '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/human_atlas/plot/referenceVec/umap/'
# d.integrated <- readRDS(datafn)
# u = d.integrated@reductions$umap@cell.embeddings
# str(u)
# meta = d.integrated@meta.data
# meta[is.na(meta$cell),'cell'] <- row.names(meta)[is.na(meta$cell)]
# u = u[meta$cell,]
u = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/human_atlas/result/referenceVec/umap/umap.rds')
meta = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/human_atlas/result/referenceVec/meta.rds')
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(gridExtra)
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

## mark GBM 'celltype' as 'GBM'
meta.bak = meta
meta[grepl('GBM',meta$study),'celltype'] <- meta[grepl('GBM',meta$study),'location'] <- meta[grepl('GBM',meta$study),'time'] <- 'GBM'
png(paste0(plotdir,'reference_umap.png'),height=800,width=2000)
p1 <- ggplot() + geom_point(data=data.frame(umap1=u[,1],umap2=u[,2],ct=meta$celltype),aes(x=umap1,y=umap2,col=ct),alpha=0.5,size=0.2) + 
  theme_classic() + xlab('UMAP1')+ylab('UMAP2')+
  theme(legend.position = 'right',legend.title = element_blank()) + 
  guides(color = guide_legend(override.aes = list(size = 4,alpha=1)))+
  theme(legend.spacing.x = unit(-0.1, 'cm'),legend.spacing.y = unit(-0.3, 'cm'))+
  scale_color_manual(values=getPalette(length(unique(meta$celltype))+1)[1:length(unique(meta$celltype))])


p2 <- ggplot() + geom_point(data=data.frame(umap1=u[,1],umap2=u[,2],location=meta$location),aes(x=umap1,y=umap2,col=location),alpha=0.5,size=0.2) + 
  theme_classic() + xlab('UMAP1')+ylab('UMAP2')+
  theme(legend.position = 'right',legend.title = element_blank()) + 
  guides(color = guide_legend(override.aes = list(size = 4,alpha=1)))+
  theme(legend.spacing.x = unit(-0.1, 'cm'),legend.spacing.y = unit(-0.3, 'cm'))+
  scale_color_manual(values=getPalette(length(unique(meta$location))+1)[1:length(unique(meta$location))])


p3 <- ggplot() + geom_point(data=data.frame(umap1=u[,1],umap2=u[,2],study=meta$study),aes(x=umap1,y=umap2,col=study),alpha=0.5,size=0.2) + 
  theme_classic() + xlab('UMAP1')+ylab('UMAP2')+
  theme(legend.position = 'right',legend.title = element_blank()) + 
  guides(color = guide_legend(override.aes = list(size = 4,alpha=1)))+
  theme(legend.spacing.x = unit(-0.1, 'cm'),legend.spacing.y = unit(-0.3, 'cm'))+
  scale_color_manual(values=getPalette(length(unique(meta$study))+1)[1:length(unique(meta$study))])

p4 <- ggplot() + geom_point(data=data.frame(umap1=u[,1],umap2=u[,2],time=meta$time),aes(x=umap1,y=umap2,col=time),alpha=0.5,size=0.2) + 
  theme_classic() + xlab('UMAP1')+ylab('UMAP2')+
  theme(legend.position = 'right',legend.title = element_blank()) + 
  guides(color = guide_legend(override.aes = list(size = 4,alpha=1)))+
  theme(legend.spacing.x = unit(-0.1, 'cm'),legend.spacing.y = unit(-0.3, 'cm'))+
  scale_color_manual(values=getPalette(length(unique(meta$time))))
grid.arrange(p1,p2,p3,p4,nrow=2)
dev.off()


### facet wrap
meta <- meta.bak
meta[grepl('GBM',meta$study),'celltype_super']  <- 'GBM'
meta[grepl('GBM',meta$study),'celltype'] <- meta[grepl('GBM',meta$study),'location'] <- meta[grepl('GBM',meta$study),'time'] <-  meta[grepl('GBM',meta$study),'study']
png(paste0(plotdir,'reference_umap_celltype.png'),height=800,width=1500)
ggplot() + geom_point(data=data.frame(umap1=u[,1],umap2=u[,2],ct=meta$celltype,ct_super=meta$celltype_super),aes(x=umap1,y=umap2,col=ct),alpha=0.5,size=0.2)+
  theme_classic() + xlab('UMAP1')+ylab('UMAP2')+
  theme(legend.position = 'right',legend.title = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 4,alpha=1)))+
  theme(legend.spacing.x = unit(-0.1, 'cm'),legend.spacing.y = unit(-0.3, 'cm'))+
  facet_wrap(~ct_super) +
  scale_color_manual(values=getPalette(length(unique(meta$celltype))))
dev.off()

###
png(paste0(plotdir,'reference_umap_celltype2.png'),height=800,width=1500)
meta <- meta.bak
meta[grepl('GBM',meta$study),'time_super']  <- meta[grepl('GBM',meta$study),'location_super'] <- meta[grepl('GBM',meta$study),'celltype_super'] <- meta[grepl('GBM',meta$study),'celltype'] <- meta[grepl('GBM',meta$study),'location'] <- meta[grepl('GBM',meta$study),'time'] <- 'GBM'
ggplot() + geom_point(data=data.frame(umap1=u[,1],umap2=u[,2],ct=meta$celltype,ct_super=meta$celltype_super),aes(x=umap1,y=umap2,col=ct),alpha=0.5,size=0.2)+
  theme_classic() + xlab('UMAP1')+ylab('UMAP2')+
  theme(legend.position = 'right',legend.title = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 4,alpha=1)))+
  theme(legend.spacing.x = unit(-0.1, 'cm'),legend.spacing.y = unit(-0.3, 'cm'))+
  facet_wrap(~ct_super) +
  scale_color_manual(values=getPalette(length(unique(meta$celltype))))
dev.off()

###
# pdf(paste0(plotdir,'reference_umap_celltype.pdf'),height=8,width=12)
# ggplot() + geom_point(data=data.frame(umap1=u[,1],umap2=u[,2],ct=meta$celltype,ct_super=meta$celltype_super),aes(x=umap1,y=umap2,col=ct),alpha=0.5,size=0.2)+
#   theme_classic() + xlab('UMAP1')+ylab('UMAP2')+
#   theme(legend.position = 'right',legend.title = element_blank()) +
#   guides(color = guide_legend(override.aes = list(size = 4,alpha=1)))+
#   theme(legend.spacing.x = unit(-0.1, 'cm'),legend.spacing.y = unit(-0.3, 'cm'))+
#   facet_wrap(~ct_super) +
#   scale_color_manual(values=getPalette(length(unique(meta$celltype))))
# dev.off()


png(paste0(plotdir,'reference_umap_location.png'),height=700,width=1050)
ggplot() + geom_point(data=data.frame(umap1=u[,1],umap2=u[,2],location=meta$location,location_super=meta$location_super),aes(x=umap1,y=umap2,col=location),alpha=0.5,size=0.2)+
  theme_classic() + xlab('UMAP1')+ylab('UMAP2')+
  theme(legend.position = 'right',legend.title = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 4,alpha=1)))+
  theme(legend.spacing.x = unit(-0.1, 'cm'),legend.spacing.y = unit(-0.3, 'cm'))+
  facet_wrap(~location_super) +
  scale_color_manual(values=getPalette(length(unique(meta$location))))
dev.off()


png(paste0(plotdir,'reference_umap_time.png'),height=700,width=1050)
ggplot() + geom_point(data=data.frame(umap1=u[,1],umap2=u[,2],time=meta$time,time_super=meta$time_super),aes(x=umap1,y=umap2,col=time),alpha=0.5,size=0.2)+
  theme_classic() + xlab('UMAP1')+ylab('UMAP2')+
  theme(legend.position = 'right',legend.title = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 4,alpha=1)))+
  theme(legend.spacing.x = unit(-0.1, 'cm'),legend.spacing.y = unit(-0.3, 'cm'))+
  facet_wrap(~time_super)+
  scale_color_manual(values=getPalette(length(unique(meta$time))))
dev.off()


png(paste0(plotdir,'reference_umap_study.png'),height=700,width=1050)
ggplot() + geom_point(data=data.frame(umap1=u[,1],umap2=u[,2],study=meta$study),aes(x=umap1,y=umap2,col=study),alpha=0.5,size=0.2)+
  theme_classic() + xlab('UMAP1')+ylab('UMAP2')+
  theme(legend.position = 'right',legend.title = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 4,alpha=1)))+
  theme(legend.spacing.x = unit(-0.1, 'cm'),legend.spacing.y = unit(-0.3, 'cm'))+
  facet_wrap(~study)+
  scale_color_manual(values=getPalette(length(unique(meta$study))+1)[1:length(unique(meta$study))])
dev.off()


### more on cell types
meta <- meta.bak
meta[grepl('GBM',meta$study),'celltype_super']  <- 'GBM'
meta[grepl('GBM',meta$study),'celltype'] <- meta[grepl('GBM',meta$study),'location'] <- meta[grepl('GBM',meta$study),'time'] <-  meta[grepl('GBM',meta$study),'study']

df = data.frame(umap1=u[,1],umap2=u[,2],ct=meta$celltype,ct_super=meta$celltype_super)

png(paste0(plotdir,'reference_umap_celltype_neuronal.png'),height=700,width=1200)
pd = df[grep('neuronal',df$ct_super),]
ggplot() + geom_point(data=pd,aes(x=umap1,y=umap2,col=ct),alpha=0.5,size=0.2)+
  theme_classic() + xlab('UMAP1')+ylab('UMAP2')+
  theme(legend.position = 'right',legend.title = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 4,alpha=1)))+
  theme(legend.spacing.x = unit(-0.1, 'cm'),legend.spacing.y = unit(-0.3, 'cm'))+
  scale_color_manual(values=getPalette(length(unique(pd$ct))))
dev.off()

pdf(paste0(plotdir,'reference_umap_celltype_neuronal.pdf'),height=8.5,width=13)
pd = df[grep('neuronal',df$ct_super),]
n = sub('\\(.*','',pd$ct)
pd = pd[order(n),]
ggplot() + geom_point(data=pd,aes(x=umap1,y=umap2,col=ct),alpha=0.5,size=0.2)+
  theme_classic() + xlab('UMAP1')+ylab('UMAP2')+
  theme(legend.position = 'right',legend.title = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 4,alpha=1)))+
  theme(legend.spacing.x = unit(-0.1, 'cm'),legend.spacing.y = unit(-0.3, 'cm'))+
  scale_color_manual(values=getPalette(length(unique(pd$ct))))
dev.off()


png(paste0(plotdir,'reference_umap_celltype_neuronal_details.png'),height=700,width=1200)
pd = df[grep('neuronal',df$ct_super),]
ggplot() + geom_point(data=pd,aes(x=umap1,y=umap2,col=ct),alpha=0.5,size=0.2)+
  theme_classic() + xlab('UMAP1')+ylab('UMAP2')+
  theme(legend.position = 'right',legend.title = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 4,alpha=1)))+
  theme(legend.spacing.x = unit(-0.1, 'cm'),legend.spacing.y = unit(-0.3, 'cm'))+
  facet_wrap(~ct)+
  scale_color_manual(values=getPalette(length(unique(pd$ct))))
dev.off()

pdf(paste0(plotdir,'reference_umap_celltype_neuronal_details.pdf'),height=7,width=12)
pd = df[grep('neuronal',df$ct_super),]
ggplot() + geom_point(data=pd,aes(x=umap1,y=umap2,col=ct),alpha=0.5,size=0.2)+
  theme_classic() + xlab('UMAP1')+ylab('UMAP2')+
  theme(legend.position = 'right',legend.title = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 4,alpha=1)))+
  theme(legend.spacing.x = unit(-0.1, 'cm'),legend.spacing.y = unit(-0.3, 'cm'))+
  facet_wrap(~ct)+
  scale_color_manual(values=getPalette(length(unique(pd$ct))))
dev.off()

png(paste0(plotdir,'reference_umap_celltype_nonNeuronal.png'),height=700,width=1200)
pd = df[!grepl('neuronal',df$ct_super),]
ggplot() + geom_point(data=pd,aes(x=umap1,y=umap2,col=ct),alpha=0.5,size=0.2)+
  theme_classic() + xlab('UMAP1')+ylab('UMAP2')+
  theme(legend.position = 'right',legend.title = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 4,alpha=1)))+
  theme(legend.spacing.x = unit(-0.1, 'cm'),legend.spacing.y = unit(-0.3, 'cm'))+
  scale_color_manual(values=getPalette(length(unique(pd$ct))))
dev.off()

pdf(paste0(plotdir,'reference_umap_celltype_nonNeuronal.pdf'),height=7,width=12)
pd = df[!grepl('neuronal',df$ct_super),]
ggplot() + geom_point(data=pd,aes(x=umap1,y=umap2,col=ct),alpha=0.5,size=0.2)+
  theme_classic() + xlab('UMAP1')+ylab('UMAP2')+
  theme(legend.position = 'right',legend.title = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 4,alpha=1)))+
  theme(legend.spacing.x = unit(-0.1, 'cm'),legend.spacing.y = unit(-0.3, 'cm'))+
  scale_color_manual(values=getPalette(length(unique(pd$ct))))
dev.off()
