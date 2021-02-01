setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/')
pdir <- './GBMonly/integrate/harmony/seuratGene/plot/'
u <- readRDS('./GBMonly/integrate/harmony/seuratGene/res/umap_embeddings.rds')
meta <- readRDS('./humanGBM/data/meta_allcell.rds')
mat <- readRDS('./humanGBM/data/combine_mat.rds')
mat <- mat[, rownames(u)]  ###
rownames(meta) <- meta$cell  ###
meta <- meta[rownames(u), ] ###


for (i in c("study",'celltype_super','celltype','time_super','time')){
  meta[is.na(meta[,i]),i] = 'unknown'
  tab = table(meta[,i])
  meta[,i] = paste0(meta[,i],'(',tab[match(meta[,i],names(tab))],')')    
}

dir.create(pdir,showWarnings = F, recursive = T)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(gridExtra)
library(scattermore)
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

for (i in colnames(meta)){
  if (!i %in% c('cell', 'cell.shortname', 'species', 'Mutation.Rate')){
    print(i)
    id = !grepl('unknown',meta[,i])  & !is.na(meta[,i])
    if (sum(id) > 0){
      p <- ggplot() + geom_scattermore(data=data.frame(umap1=u[id,1],umap2=u[id,2],feature=meta[id,i]),aes(x=umap1,y=umap2,col=feature),alpha=0.2,size=0.1) + 
        theme_classic() + xlab('UMAP1')+ylab('UMAP2')+
        theme(legend.position = 'right',legend.title = element_blank()) + 
        guides(color = guide_legend(override.aes = list(size = 4,alpha=1)))+
        theme(legend.spacing.x = unit(-0.1, 'cm'),legend.spacing.y = unit(-0.3, 'cm'))+
         scale_color_manual(values=getPalette(length(unique(meta[,i]))+1)[1:length(unique(meta[,i]))])
      ggsave(paste0(pdir,'umap_', i, '.png'),p,height=3.5 + 0.5*(length(unique(meta[,i]))/8),width=6+length(unique(meta[,i]))/12,dpi=200)
    }
  }
}
  

id.gbm <- grepl('GBM', meta$study) 
for (i in colnames(meta)){
  if (!i %in% c('cell', 'cell.shortname', 'species', 'Mutation.Rate')){
    print(i)
    id = !grepl('unknown',meta[,i])  & !is.na(meta[,i])
    id <- id & id.gbm
    if (sum(id) > 0){
      p <- ggplot() + geom_scattermore(data=data.frame(umap1=u[id,1],umap2=u[id,2],feature=meta[id,i]),aes(x=umap1,y=umap2,col=feature),alpha=0.2,size=0.1) + 
        theme_classic() + xlab('UMAP1')+ylab('UMAP2')+
        theme(legend.position = 'right',legend.title = element_blank()) + 
        guides(color = guide_legend(override.aes = list(size = 4,alpha=1)))+
        theme(legend.spacing.x = unit(-0.1, 'cm'),legend.spacing.y = unit(-0.3, 'cm'))+
         scale_color_manual(values=getPalette(length(unique(meta[,i]))+1)[1:length(unique(meta[,i]))])
      ggsave(paste0(pdir,'umap_GBM_', i, '.png'),p,height=3.5 + 0.5*(length(unique(meta[,i]))/12),width=6+length(unique(meta[,i]))/12,dpi=200)
    }
  }
}

for (i in colnames(meta)){
  if (!i %in% c('cell', 'cell.shortname', 'species', 'Mutation.Rate')){
    print(i)
    id = !grepl('unknown',meta[,i])  & !is.na(meta[,i])
    id <- id & (!id.gbm)
    if (sum(id) > 0){
      p <- ggplot() + geom_scattermore(data=data.frame(umap1=u[id,1],umap2=u[id,2],feature=meta[id,i]),aes(x=umap1,y=umap2,col=feature),alpha=0.2,size=0.1) + 
        theme_classic() + xlab('UMAP1')+ylab('UMAP2')+
        theme(legend.position = 'right',legend.title = element_blank()) + 
        guides(color = guide_legend(override.aes = list(size = 4,alpha=1)))+
        theme(legend.spacing.x = unit(-0.1, 'cm'),legend.spacing.y = unit(-0.3, 'cm'))+
         scale_color_manual(values=getPalette(length(unique(meta[,i]))+1)[1:length(unique(meta[,i]))])
      ggsave(paste0(pdir,'umap_atlas_', i, '.png'),p,height=3.5 + 0.5*(length(unique(meta[,i]))/12),width=6+length(unique(meta[,i]))/12,dpi=200)
    }
  }
}
   
   
## farcet study, celltype
for (i in  c('study','celltype')){
  id = !grepl('unknown',meta[,i])  & !is.na(meta[,i])
    id <- id & (!id.gbm)
    if (sum(id) > 0){
      p <- ggplot() + geom_scattermore(data=data.frame(umap1=u[id,1],umap2=u[id,2],feature=meta[id,i]),aes(x=umap1,y=umap2,col=feature),alpha=0.8,size=10) + 
        theme_classic() + xlab('UMAP1')+ylab('UMAP2')+
        theme(legend.position = 'right',legend.title = element_blank()) + 
        guides(color = guide_legend(override.aes = list(size = 4,alpha=1)))+
        theme(legend.spacing.x = unit(-0.1, 'cm'),legend.spacing.y = unit(-0.3, 'cm'))+
        facet_wrap(~feature, nrow = round(sqrt(length(unique(meta[id,i]))))) +
        theme(legend.position = 'none')
      if (length(unique(meta[id,i]))>8){
        p <- p + scale_color_manual(values=getPalette(length(unique(meta[,i]))+1)[1:length(unique(meta[,i]))])
      } else {
        p <- p + scale_color_brewer(palette = 'Set1')
      }
      ggsave(paste0(pdir,'umap_atlas_', i, '_facet.png'),p,height=3.5 + (length(unique(meta[,i]))/6),width=6+length(unique(meta[,i]))/3,dpi=200)
    }
}




  p <- ggplot() + geom_scattermore(data=data.frame(umap1=u[,1],umap2=u[,2],feature=ifelse(id.gbm, 'GBM', 'atlas')),aes(x=umap1,y=umap2,col=feature),alpha=0.2,size=0.1) + 
    theme_classic() + xlab('UMAP1')+ylab('UMAP2')+
    theme(legend.position = 'right',legend.title = element_blank()) + 
    guides(color = guide_legend(override.aes = list(size = 4,alpha=1)))+
    theme(legend.spacing.x = unit(-0.1, 'cm'),legend.spacing.y = unit(-0.3, 'cm'))+
     scale_color_manual(values=c('orange', 'steelblue'))
  ggsave(paste0(pdir,'umap_GBM_atlas.png'),p,height=4, width=5,dpi=200)



### plot marker genes
plotfunc <- function(gene){
  v <-  mat[gene,]
  v[v > quantile(v, 0.98)] <- quantile(v, 0.98)
  v[v < quantile(v, 0.02)] <- quantile(v, 0.02)  
  p <- ggplot() + geom_scattermore(data=data.frame(umap1=u[,1],umap2=u[,2],expr = v),aes(x=umap1,y=umap2,col=expr),alpha=1,size=0.2) + 
  theme_classic() + xlab('UMAP1')+ylab('UMAP2')+
  ggtitle(gene)+ 
  scale_color_gradientn(colors=brewer.pal(11,'Spectral')) 
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

ggsave(paste0(pdir,'umap_marker_gene.png'),grid.arrange(grobs=plist,nrow=4),
       height=10,width=14,
       dpi=200)

plotfunc2 <- function(gene,ct){
  v <-  mat[gene,]
  v[v > quantile(v, 0.98)] <- quantile(v, 0.98)
  v[v < quantile(v, 0.02)] <- quantile(v, 0.02)  
  p <- ggplot() + geom_scattermore(data=data.frame(umap1=u[,1],umap2=u[,2],expr = v),aes(x=umap1,y=umap2,col=expr),alpha=1,size=0.2) + 
  theme_classic() + xlab('UMAP1')+ylab('UMAP2')+
  ggtitle(paste0(ct,':',gene))+ 
  scale_color_gradientn(colors=brewer.pal(11,'Spectral')) 
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
  ggsave(paste0(pdir,paste0('umap_marker_gene_',colnames(tb)[i],'.png')),grid.arrange(grobs=plist,nrow=4,ncol=4),
       height=10,width=14,
       dpi=200)
}

