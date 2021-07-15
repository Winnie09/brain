rdir = '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/humanatlascortexGBM/seurat/res/'
res = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/humanatlascortexGBM/seurat/res/project.rds')
loc = read.csv('/home-4/whou10@jhu.edu/work-zfs/whou10/GBM/doc/metas_20200329.csv', as.is = T)
loc <- loc[loc$Pathology=='GBM' & loc$Tumor.Grade=='IV' & loc$Treatment=='Untreated',]

## atlas cell cluster belongs to which patient
tab <- table(res@meta.data$predicted.celltype,sub('_.*','',rownames(res@meta.data)))
write.csv(tab, '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/humanatlascortexGBM/seurat/res/GBM_predicted_cellnumber.csv')

tab <- t(tab)/colSums(tab)
write.csv(tab, '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/humanatlascortexGBM/seurat/res/GBM_predicted_cellproportion.csv')

clu <- table(res@meta.data$predicted.celltype)
selclu <- names(clu)[clu>10]

tab <- tab[intersect(rownames(tab),loc[,1]),selclu]
write.csv(tab, paste0(rdir, 'GBM_predicted_cellproportion_cluster10morecells.csv'))

## location of each sample
sploc <- loc[match(rownames(tab),loc[,1]),'Location']

## test one location vs the other location: greater
pval <- sapply(unique(sploc),function(i) {
  sapply(colnames(tab),function(j) {
    t.test(tab[sploc==i,j],tab[sploc!=i,j],alternative = 'greater')$p.value
  })
})
write.csv(pval, paste0(rdir, 'pval_loc_greater.csv'))

## test one location vs the other location: less
pval <- sapply(unique(sploc),function(i) {
  sapply(colnames(tab),function(j) {
    t.test(tab[sploc==i,j],tab[sploc!=i,j],alternative = 'less')$p.value
  })
})
write.csv(pval, paste0(rdir, 'pval_loc_less.csv'))
which(pval < 0.05,arr.ind=T)

tmp <- tab[,'13']
names(tmp) <- sploc
tmp

tmp <- tab[,'10']
names(tmp) <- sploc
tmp

tab2 = cbind(sploc, tab)
write.csv(tab2, paste0(rdir, 'GBM_predicted_cellproportion_cluster10morecells_with_loc.csv'))
