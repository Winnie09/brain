library(Seurat)
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin')
d = readRDS('./brain/atlas/result/atlasOnly_homolog_sub/integrated.rds')
# fullmeta = readRDS('./brain/atlas/old/20200212/result/full/meta.rds')
mat = as.matrix(d@assays$RNA@counts)
s = gsub(':.*','', colnames(mat))

u = readRDS('./brain/atlas/result/atlasOnly_homolog_sub/umap/umap.rds')

meta = data.frame(cell = rownames(u), study = sub(':.*','',rownames(u)))
meta = cbind(meta, 
             celltype_super = fullmeta[match(meta$cell, fullmeta$cell),'celltype_super'],
             species = fullmeta[match(meta$cell, fullmeta$cell),'species'],
             time_super = fullmeta[match(meta$cell, fullmeta$cell),'time_super'])

meta[grepl('2017_Nowakowski_Science',meta$study),'species'] = 'human'
meta[grepl('2016_Tasic_NatNeuro',meta$study),'species'] = 'mouse'
meta[grepl('2015_Zeisel_Science',meta$study),'species'] = 'mouse'
meta[grepl('2017_Nowakowski_Science',meta$study),'time_super'] = 'developing'
meta[grepl('2016_Tasic_NatNeuro',meta$study),'time_super'] = 'adult'
meta[grepl('2015_Zeisel_Science',meta$study),'time_super'] = 'adult'
meta[is.na(meta$species),'species'] <- 'mouse'


for (i in c("study",'celltype_super','species','time_super')){
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
p3 <- ggplot() + geom_point(data=data.frame(umap1=u[id,1],umap2=u[id,2],study=meta$species[id]),aes(x=umap1,y=umap2,col=study),alpha=0.2,size=0.1) + 
  theme_classic() + xlab('UMAP1')+ylab('UMAP2')+
  theme(legend.position = 'right',legend.title = element_blank()) + 
  guides(color = guide_legend(override.aes = list(size = 4,alpha=1)))+
  theme(legend.spacing.x = unit(-0.1, 'cm'),legend.spacing.y = unit(-0.3, 'cm'))+
   scale_color_manual(values=getPalette(length(unique(meta$species))+1)[1:length(unique(meta$species))])

id2 = !grepl('unknown',meta$time_super)
p4 <- ggplot() + geom_point(data=data.frame(umap1=u[id2,1],umap2=u[id2,2],study=meta$time_super[id2]),aes(x=umap1,y=umap2,col=study),alpha=0.2,size=0.1) + 
  theme_classic() + xlab('UMAP1')+ylab('UMAP2')+
  theme(legend.position = 'right',legend.title = element_blank()) + 
  guides(color = guide_legend(override.aes = list(size = 4,alpha=1)))+
  theme(legend.spacing.x = unit(-0.1, 'cm'),legend.spacing.y = unit(-0.3, 'cm'))+
   scale_color_manual(values=getPalette(length(unique(meta$time_super))+1)[1:length(unique(meta$time_super))])

pdf(paste0(plotdir,'umap.pdf'),height=6,width=15)
grid.arrange(p1,p2,p3,p4,nrow=2, layout_matrix = rbind(c(1,1,2,2,2),c(3,3,4,4,5)))
dev.off()


png(paste0(plotdir,'umap_time.png'),height=300,width=500)
ggplot() + geom_point(data=data.frame(umap1=u[id2,1],umap2=u[id2,2],study=meta$time_super[id2]),aes(x=umap1,y=umap2,col=study),alpha=0.2,size=0.1) + 
  theme_classic() + xlab('UMAP1')+ylab('UMAP2')+
  theme(legend.position = 'right',legend.title = element_blank()) + 
  guides(color = guide_legend(override.aes = list(size = 4,alpha=1)))+
  theme(legend.spacing.x = unit(-0.1, 'cm'),legend.spacing.y = unit(-0.3, 'cm'))+
   scale_color_manual(values=getPalette(length(unique(meta$time_super))+1)[1:length(unique(meta$time_super))])+
  facet_wrap(~study)
dev.off()


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
         'POU5F1', 'SOX2', 'KLF4', 'NANOG','MYC') ## embryonic
gvec = intersect(gvec, rownames(mat))

plist = list()
for (g in gvec){
  plist[[g]] = plotfunc(g)
}
pdf(paste0(plotdir,'umap_marker_gene.pdf'),height=7,width=12)
grid.arrange(grobs=plist,nrow=4)
dev.off()
png(paste0(plotdir,'umap_marker_gene.png'),height=700,width=1200)
grid.arrange(grobs=plist,nrow=4)
dev.off()



