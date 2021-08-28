library(Seurat)
rdir = '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/humanatlascortexGBM/seurat/res/'
pdir = '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/humanatlascortexGBM/seurat/plot/'

# res <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/humanatlascortexGBM/seurat/res/project.rds')

## get GBM projected predicted celltype and the annotation
# gbm.pred.ct <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/humanatlascortexGBM/seurat/res/gbm_predict_celltpye_anno.rds')




res <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/humanatlascortexGBM/seurat/res/project_with_celltype_anno.rds')


normalclu <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_myelymreference/res/normalclu.rds')

clu <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/integrate/seurat/36nonNormalSeuratGene2000/res/cluster.rds')

cell <- names(clu)[!clu %in% normalclu]

meta <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/integrate/seurat/36nonNormalSeuratGene2000/res/meta.rds')

meta <- meta[,c('Sample.ID','Location')]
loc <- unique(meta)
# loc = read.csv('/home-4/whou10@jhu.edu/work-zfs/whou10/GBM/doc/metas_20200329.csv', as.is = T)
# loc <- loc[loc$Pathology=='GBM' & loc$Tumor.Grade=='IV' & loc$Treatment=='Untreated',]

## atlas cell cluster belongs to which patient

resmeta <- res@meta.data
resmeta <- resmeta[cell,]

## predicted cluster by patient
tab <- table(resmeta$predicted.celltype,sub('_.*','',rownames(resmeta)))
write.csv(tab, '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/humanatlascortexGBM/seurat/res/GBM_predicted_cluster_by_cellnumber_table.csv')

tab <- t(tab)/colSums(tab)
write.csv(tab, '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/humanatlascortexGBM/seurat/res/GBM_predicted_cluster_by_cellprop_table.csv')

clu <- table(resmeta$predicted.celltype)
selclu <- names(clu)[clu>10]
tab <- tab[intersect(rownames(tab),loc[,1]),selclu]
write.csv(tab, paste0(rdir, 'GBM_predicted_cluster_by_cellproportion_cluster10morecells.csv'))

## predicted cell type by patient
tab <- table(resmeta$predicted.celltype_anno,sub('_.*','',rownames(resmeta)))
write.csv(tab, '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/humanatlascortexGBM/seurat/res/GBM_predicted_celltype_by_cellnumber_table.csv')

tab <- t(tab)/colSums(tab)
write.csv(tab, '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/humanatlascortexGBM/seurat/res/GBM_predicted_celltype_by_cellprop_table.csv')

ct <- table(resmeta$predicted.celltype_anno)
selct <- names(ct)[ct>10]
tab <- tab[intersect(rownames(tab),loc[,1]),selct]
write.csv(tab, paste0(rdir, 'GBM_predicted_celltype_by_cellprop_cluster10morecells.csv'))


## location of each sample

pd = t(tab)
sploc <- loc[match(colnames(pd),loc[,1]),'Location']
colann = data.frame(location = sploc, stringsAsFactors = F)
rownames(colann) = colnames(pd)
  
library(pheatmap)
pdf(paste0(pdir, 'GBM_predictedd_celltype_by_cellprop_heatmap.pdf'), height = 3, width = 8)
pheatmap(pd, annotation_col = colann)
dev.off()

## save a pvalue table with samples marked locations
tab2 = cbind(sploc, tab)
write.csv(tab2, paste0(rdir, 'GBM_predicted_celltype_cellproportion_cluster10morecells_with_loc.csv'))

### ====== ###
### test   ###
### ====== ###
## test one location vs the other location: greater
pval <- sapply(unique(sploc),function(i) {
  sapply(colnames(tab),function(j) {
    t.test(tab[sploc==i,j],tab[sploc!=i,j],alternative = 'greater')$p.value
  })
})
write.csv(pval, paste0(rdir, 'pval_loc_greater.csv'))
saveRDS(pval, paste0(rdir, 'pval_loc_greater.rds'))


library(RColorBrewer)
library(reshape2)
pd <- melt(pval)
pd$sig <- ifelse(pd$value < 0.05,'*','')
pd$value <- -log10(pd$value)
pdf(paste0(pdir, 'GBM_predictedd_celltype_by_location_pvalue_greater_heatmap.pdf'), height = 3.5, width = 3.5)
ggplot(pd,aes(x=Var1,y=Var2,fill=value,label=sig)) + geom_tile() + geom_text() + theme_classic() + scale_fill_gradient2(low='blue',mid='white',high='red',midpoint=-log10(0.05)) + theme(legend.position = 'bottom') + xlab('') + ylab('') + labs(fill="-log10 pvalue") + coord_flip()+
  theme(axis.text.x = element_text(angle=45, colour = 'black', hjust = 1), 
        axis.text.y = element_text(color='black'))
dev.off()




## test one location vs the other location: less
pval <- sapply(unique(sploc),function(i) {
  sapply(colnames(tab),function(j) {
    t.test(tab[sploc==i,j],tab[sploc!=i,j],alternative = 'less')$p.value
  })
})
write.csv(pval, paste0(rdir, 'pval_loc_less.csv'))
saveRDS(pval, paste0(rdir, 'pval_loc_less.rds'))

pd <- melt(pval)
pd$sig <- ifelse(pd$value < 0.05,'*','')
pd$value <- -log10(pd$value)
pdf(paste0(pdir, 'GBM_predictedd_celltype_by_location_pvalue_less_heatmap.pdf'), height = 3.5, width = 3.5)
ggplot(pd,aes(x=Var1,y=Var2,fill=value,label=sig)) + geom_tile() + geom_text() + theme_classic() + scale_fill_gradient2(low='blue',mid='white',high='red',midpoint=-log10(0.05)) + theme(legend.position = 'bottom') + xlab('') + ylab('') + labs(fill="-log10 pvalue") + coord_flip()+
  theme(axis.text.x = element_text(angle=45, colour = 'black', hjust = 1), 
        axis.text.y = element_text(color='black'))
dev.off()


which(pval < 0.05,arr.ind=T)


