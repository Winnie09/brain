library(here)
setwd(here())
mat =  readRDS('/home-4/whou10@jhu.edu/work-zfs/whou10/GBM/data/proc/tumor/all/matrix/count.rds')
str(mat)
allp <- sub('_.*', '', colnames(mat))
table(allp)
id <- which(!allp %in% c('GBM076', 'GBM077', 'GBM086', 'GBM089', 'GBM090', 'GBM094', 'GBM047')) ## six normal and a not-GBM related patients
mat <- mat[,id]
str(mat)
allp <- allp[id]
table(allp)
saveRDS(mat, 'atlasGBM/GBMonly/data/36nonNormal_combined_count_mat.rds')
mat <- log2(t(t(mat)/(colSums(mat) / 1e6) + 1))
summary(colSums(mat))
saveRDS(mat, 'atlasGBM/GBMonly/data/36nonNormal_combined_log2norm_mat.rds')

dlist = list()
for (i in unique(allp)){
  dlist[[i]] <- mat[, allp == i]
}
saveRDS(dlist, 'atlasGBM/GBMonly/data/36nonNormal_combined_log2norm_dlist.rds')


meta <- read.csv('/home-4/whou10@jhu.edu/work-zfs/whou10/GBM/doc/metas_20200329.csv',as.is=T)
cm <- data.frame(cell=colnames(mat),study='GBM',meta[match(sub('_.*','',colnames(mat)),meta[,1]),])
saveRDS(cm, 'atlasGBM/GBMonly/data/36nonNormal_combined_meta.rds')


### 6 normal samples
cnt = readRDS('/home-4/whou10@jhu.edu/work-zfs/whou10/GBM/data/proc/tumor/matrix/count.rds')
meta.cnt = read.csv('/home-4/whou10@jhu.edu/work-zfs/whou10/GBM/doc/metas_20200329.csv', as.is = T)
s = meta.cnt[meta.cnt$Pathology == 'Normal', 1]
as = sub('_.*', '', colnames(cnt))
cnt = cnt[, as %in% s]

mat <- log2(t(t(cnt)/(colSums(cnt) / 1e6) + 1))

as = sub('_.*', '', colnames(mat))

dlist = list()
for (i in unique(as)){
  dlist[[i]] <- mat[, as == i]
}
saveRDS(dlist, 'atlasGBM/GBMonly/data/6Normal_combined_log2norm_dlist.rds')


