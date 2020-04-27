library(Seurat)
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin')
d = readRDS('./brain/atlas/result/atlasOnly_homolog_sub/integrated.rds')
fullmeta = readRDS('./brain/atlas/old/20200212/result/full/meta.rds')
mat = as.matrix(d@assays$RNA@counts)
s = gsub(':.*','', colnames(mat))

u = readRDS('./brain/atlas/result/atlasOnly_homolog_sub/umap/umap.rds')

meta = data.frame(cell = rownames(u), study = sub(':.*','',rownames(u)))
meta = cbind(meta, 
             celltype_super = fullmeta[match(meta$cell, fullmeta$cell),'celltype_super'],
             species = fullmeta[match(meta$cell, fullmeta$cell),'species'],
             time_super = fullmeta[match(meta$cell, fullmeta$cell),'time_super'],
             celltype = fullmeta[match(meta$cell, fullmeta$cell),'celltype'],
             time = fullmeta[match(meta$cell, fullmeta$cell),'time'])


meta$time <- sapply(1:nrow(meta), function(i){
  if (grepl('unknown',meta[i,'time'])){
    meta[i,'time_super']
  } else {
    meta[i,'time']
  }
})


meta[grepl('2017_Nowakowski_Science',meta$study),'species'] = 'human'
meta[grepl('2016_Tasic_NatNeuro',meta$study),'species'] = 'mouse'
meta[grepl('2015_Zeisel_Science',meta$study),'species'] = 'mouse'
meta[grepl('2017_Nowakowski_Science',meta$study),'time_super'] = 'developing'
meta[grepl('2016_Tasic_NatNeuro',meta$study),'time_super'] = 'adult'
meta[grepl('2015_Zeisel_Science',meta$study),'time_super'] = 'adult'
meta[is.na(meta$species),'species'] <- 'mouse'


for (i in c("study",'celltype_super','celltype','species','time_super','time')){
  meta[is.na(meta[,i]),i] = 'unknown'
  tab = table(meta[,i])
  meta[,i] = paste0(meta[,i],'(',tab[match(meta[,i],names(tab))],')')    
  if (i == 'study')
  for (j in 1:nrow(meta)){
    if (grepl('2018_Rosenberg',meta[j,2])){
      meta[j,2] = sub(')', paste0('/','156049',')'),meta[j,2] )
    } else if (grepl('2018_Zeisel',meta[j,2])) {
      meta[j,2] = sub(')', paste0('/','492949',')'),meta[j,2] )
    } else if (grepl('allen',meta[j,2])){
      meta[j,2] = sub(')', paste0('/','47509',')'),meta[j,2] )
    }
  }
}


plotdir = './brain/atlas/plot/atlasOnly_homolog_sub/'
dir.create(plotdir,showWarnings = F, recursive = T)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(gridExtra)
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
p1 <- ggplot() + geom_point(data=data.frame(umap1=u[,1],umap2=u[,2],study=meta$study),aes(x=umap1,y=umap2,col=study),alpha=0.2,size=0.1) + 
  theme_classic() + xlab('UMAP1')+ylab('UMAP2')+
  theme(legend.position = 'right',legend.title = element_blank()) + 
  guides(color = guide_legend(override.aes = list(size = 4,alpha=1)))+
  theme(legend.spacing.x = unit(-0.1, 'cm'),legend.spacing.y = unit(-0.3, 'cm'))+
   scale_color_manual(values=getPalette(length(unique(meta$study))+1)[1:length(unique(meta$study))])
id1 = !grepl('unknown',meta$celltype_super) 

p2 <- ggplot() + geom_point(data=data.frame(umap1=u[id1,1],umap2=u[id1,2],study=meta$celltype_super[id1]),aes(x=umap1,y=umap2,col=study),alpha=0.2,size=0.1) + 
  theme_classic() + xlab('UMAP1')+ylab('UMAP2')+
  theme(legend.position = 'right',legend.title = element_blank()) + 
  guides(color = guide_legend(override.aes = list(size = 4,alpha=1)))+
  theme(legend.spacing.x = unit(-0.1, 'cm'),legend.spacing.y = unit(-0.3, 'cm'))+
   scale_color_manual(values=getPalette(length(unique(meta$celltype_super))+1)[1:length(unique(meta$celltype_super))])
p3 <- ggplot() + geom_point(data=data.frame(umap1=u[id1,1],umap2=u[id1,2],study=meta$species[id1]),aes(x=umap1,y=umap2,col=study),alpha=0.2,size=0.1) + 
  theme_classic() + xlab('UMAP1')+ylab('UMAP2')+
  theme(legend.position = 'right',legend.title = element_blank()) + 
  guides(color = guide_legend(override.aes = list(size = 4,alpha=1)))+
  theme(legend.spacing.x = unit(-0.1, 'cm'),legend.spacing.y = unit(-0.3, 'cm'))+
   scale_color_manual(values=c('orange','steelblue'))

id2 = !grepl('unknown',meta$time_super)
p4 <- ggplot() + geom_point(data=data.frame(umap1=u[id2,1],umap2=u[id2,2],study=meta$time_super[id2]),aes(x=umap1,y=umap2,col=study),alpha=0.2,size=0.1) + 
  theme_classic() + xlab('UMAP1')+ylab('UMAP2')+
  theme(legend.position = 'right',legend.title = element_blank()) + 
  guides(color = guide_legend(override.aes = list(size = 4,alpha=1)))+
  theme(legend.spacing.x = unit(-0.1, 'cm'),legend.spacing.y = unit(-0.3, 'cm'))+
   scale_color_manual(values=getPalette(length(unique(meta$time_super))+1)[1:length(unique(meta$time_super))])

pdf(paste0(plotdir,'umap.pdf'),height=7,width=17.5)
grid.arrange(p1,p2,p3,p4,nrow=2, layout_matrix = rbind(c(1,1,2,2,2),c(3,3,4,4,5)))
dev.off()

ggsave(paste0(plotdir,'umap.png'),
       grid.arrange(p1,p2,p3,p4,nrow=2, layout_matrix = rbind(c(1,1,2,2,2),c(3,3,4,4,5))),
       height=7,width=17.5,
       dpi=200)


pp <- ggplot() + geom_point(data=data.frame(umap1=u[id1,1],umap2=u[id1,2],study=meta$species[id1]),aes(x=umap1,y=umap2,col=study),alpha=0.2,size=0.1) + 
  theme_classic() + xlab('UMAP1')+ylab('UMAP2')+
  theme(legend.position = 'right',legend.title = element_blank()) + 
  guides(color = guide_legend(override.aes = list(size = 4,alpha=1)))+
  theme(legend.spacing.x = unit(-0.1, 'cm'),legend.spacing.y = unit(-0.3, 'cm'))+
   scale_color_manual(values=c('orange','steelblue')) + facet_wrap(~study)
ggsave(paste0(plotdir,'umap_species.png'),pp,height=4,width=10,dpi=200)

p <- ggplot() + geom_point(data=data.frame(umap1=u[id2,1],umap2=u[id2,2],study=meta$time_super[id2]),aes(x=umap1,y=umap2,col=study),alpha=0.2,size=0.1) + 
  theme_classic() + xlab('UMAP1')+ylab('UMAP2')+
  theme(legend.position = 'right',legend.title = element_blank()) + 
  guides(color = guide_legend(override.aes = list(size = 4,alpha=1)))+
  theme(legend.spacing.x = unit(-0.1, 'cm'),legend.spacing.y = unit(-0.3, 'cm'))+
   scale_color_manual(values=getPalette(length(unique(meta$time_super))+1)[1:length(unique(meta$time_super))])+
  facet_wrap(~study)
ggsave(paste0(plotdir,'umap_time.png'),p,height=4,width=10,dpi=200)


p <- ggplot() + geom_point(data=data.frame(umap1=u[id1,1],umap2=u[id1,2],study=meta$celltype_super[id1]),aes(x=umap1,y=umap2,col=study),alpha=0.2,size=0.1) + 
  theme_classic() + xlab('UMAP1')+ylab('UMAP2')+
  theme(legend.position = 'right',legend.title = element_blank()) + 
  guides(color = guide_legend(override.aes = list(size = 4,alpha=1)))+
  theme(legend.spacing.x = unit(-0.1, 'cm'),legend.spacing.y = unit(-0.3, 'cm'))+
   scale_color_manual(values=getPalette(length(unique(meta$celltype_super))+1)[1:length(unique(meta$celltype_super))]) + facet_wrap(~study, ncol=6)
ggsave(paste0(plotdir,'umap_celltype.png'),p,height=8,width=18,dpi=200)

id = grep('neuronal', meta[,'celltype_super'])
p <- ggplot() + geom_point(data=data.frame(umap1=u[id,1],umap2=u[id,2],study=meta$celltype[id]),aes(x=umap1,y=umap2,col=study),alpha=0.5,size=0.1) + 
  theme_classic() + xlab('UMAP1')+ylab('UMAP2')+
  theme(legend.position = 'none',legend.title = element_blank()) + 
  guides(color = guide_legend(override.aes = list(size = 4,alpha=1)))+
  theme(legend.spacing.x = unit(-0.1, 'cm'),legend.spacing.y = unit(-0.3, 'cm'))+
   scale_color_manual(values=getPalette(length(unique(meta$celltype))+1)[1:length(unique(meta$celltype))]) + facet_wrap(~study, ncol=6)
ggsave(paste0(plotdir,'umap_neuronal.png'),p,height=8,width=18,dpi=200)

id = grep('neuronal', meta[,'celltype_super'])
p <- ggplot() + geom_point(data=data.frame(umap1=u[id,1],umap2=u[id,2],study=meta$celltype[id]),aes(x=umap1,y=umap2,col=study),alpha=0.2,size=0.1) + 
  theme_classic() + xlab('UMAP1')+ylab('UMAP2')+
  theme(legend.position = 'none',legend.title = element_blank()) + 
  guides(color = guide_legend(override.aes = list(size = 4,alpha=1)))+
  theme(legend.spacing.x = unit(-0.1, 'cm'),legend.spacing.y = unit(-0.3, 'cm'))+
   scale_color_manual(values=getPalette(length(unique(meta$celltype))+1)[1:length(unique(meta$celltype))]) + facet_wrap(~study, ncol=6)
ggsave(paste0(plotdir,'umap_neuronal.png'),p,height=8,width=17,dpi=200)
 
p <- ggplot() + geom_point(data=data.frame(umap1=u[,1],umap2=u[,2],study=meta$time),aes(x=umap1,y=umap2,col=study),alpha=0.2,size=0.1) + 
  theme_classic() + xlab('UMAP1')+ylab('UMAP2')+
  theme(legend.position = 'none',legend.title = element_blank()) + 
  guides(color = guide_legend(override.aes = list(size = 4,alpha=1)))+
  theme(legend.spacing.x = unit(-0.1, 'cm'),legend.spacing.y = unit(-0.3, 'cm'))+
   scale_color_manual(values=getPalette(length(unique(meta$celltype))+1)[1:length(unique(meta$celltype))]) + facet_wrap(~study, ncol=5)
ggsave(paste0(plotdir,'umap_time_detail.png'),p,height=9,width=12,dpi=200)


### plot marker genes
plotfunc <- function(gene){
  p <- ggplot() + geom_point(data=data.frame(umap1=u[,1],umap2=u[,2],expr = mat[gene,]),aes(x=umap1,y=umap2,col=expr),alpha=1,size=0.2) + 
  theme_classic() + xlab('UMAP1')+ylab('UMAP2')+
  ggtitle(gene)+ 
  scale_color_gradientn(colors=brewer.pal(9,'PuBu')) 
  return(p)
}

gvec = c('CD40','CD68','CD11b','CD45',  ## microglia
         'ALDH1L1','GFAP','EAA1','GLAST','EAA2','GLT1','GLT-1', ## astrocyte
         'MBP','OLIG1','OLIG2', ## oligodendrocytes
         'MAP2', 'NeuN', 'PSD95', ## mature neurons
         'SOX10', ## schwann cell precursor, myelinating schwann cell
         'SOX2','HES1', ##  neuroepithelial
         'GAT1', 'GAD65','GAD67', ## GABAergic neurons
         'TH','DAT', ## dopaminnergic neuron
         'TPH',
         'ChAT',
         'POU5F1', 'SOX2', 'KLF4', 'NANOG','MYC','OCT4','OCT') ## embryonic
gvec = intersect(gvec, rownames(mat))
# gvec = sapply(gvec, function(i) rownames(mat)[grep(i, rownames(mat))])

plist = list()
for (g in gvec){
  plist[[g]] = plotfunc(g)
}

ggsave(paste0(plotdir,'umap_marker_gene.png'),grid.arrange(grobs=plist,nrow=4),
       height=10,width=14,
       dpi=200)

plotfunc2 <- function(gene,ct){
  p <- ggplot() + geom_point(data=data.frame(umap1=u[,1],umap2=u[,2],expr = mat[gene,]),aes(x=umap1,y=umap2,col=expr),alpha=1,size=0.2) + 
  theme_classic() + xlab('UMAP1')+ylab('UMAP2')+
  ggtitle(paste0(ct,':',gene))+ 
  scale_color_gradientn(colors=brewer.pal(9,'PuBu')) 
  return(p)
}
tb = read.table('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/doc/Marker_Genes.csv',header=T,sep=',',as.is=T)
for (i in seq(1,ncol(tb))){
  print(i)
  plist = list()
  gvec = intersect(tb[,i], rownames(mat))
  for (g in gvec){
    plist[[g]] = plotfunc(g)
  }
  ggsave(paste0(plotdir,paste0('umap_marker_gene_',colnames(tb)[i],'.png')),grid.arrange(grobs=plist,nrow=4,ncol=4),
       height=10,width=14,
       dpi=200)
}
  
# 
# g = 'POU5F1'
# pd <- sapply(unique(meta$study), function(i){
#  cbind(expr = mat[g, as.character(meta$study) == i], study = i)
# })
# pd = do.call(rbind, pd)
# pd = data.frame(expr = as.numeric(pd[,1]), study = factor(pd[,2]))
# ggplot() + geom_boxplot(data = pd, aes(x=study, y = expr),  alpha = 0.3) +
#   theme(legend.position = 'none')+
#   theme(axis.text.x = element_text(angle=90)) +
#   theme_classic() + coord_flip()
# 
# table(meta[grep('2016_La', as.character(meta$study)),'celltype_super'])
