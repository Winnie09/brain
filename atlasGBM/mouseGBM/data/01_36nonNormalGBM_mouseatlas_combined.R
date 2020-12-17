library(here)
setwd(here())
gbm.meta <- readRDS('atlasGBM/GBMonly/data/36nonNormal_combined_meta.rds')
atlas.meta <- readRDS('atlasGBM/mouseatlas/data/meta_allcell.rds')
colnames(gbm.meta)[colnames(gbm.meta) == 'Age'] <- 'age'
colnames(gbm.meta)[colnames(gbm.meta) == 'Sex'] <- 'sex'
colnames(gbm.meta)[colnames(gbm.meta) == 'Grade'] <- 'grade'
colnames(atlas.meta)[colnames(atlas.meta) == 'gender'] <- 'sex'
gbm.meta$time <- 'humanAdult'
meta <- reshape::merge_recurse(list(atlas.meta, gbm.meta))
saveRDS(meta, 'atlasGBM/mouseGBM/data/36nonNormalGBM_mouseatlas_combined_meta.rds')

gbm <- readRDS('atlasGBM/GBMonly/data/36nonNormal_combined_log2norm_dlist.rds')
atlas <- readRDS('atlasGBM/mouseatlas/data/dlist.rds')
dlist <- c(atlas, gbm)
intgn <- sapply(dlist,row.names)
intgn <- table(unlist(intgn))
intgn <- names(intgn)[intgn==length(dlist)]
for (i in names(dlist)){
  dlist[[i]] <- dlist[[i]][intgn, ]
}
saveRDS(dlist, 'atlasGBM/mouseGBM/data/36nonNormalGBM_mouseatlas_combined_log2norm_dlist.rds')

mat <- do.call(cbind, dlist)
saveRDS(mat, 'atlasGBM/mouseGBM/data/36nonNormalGBM_mouseatlas_combined_log2norm_mat.rds')



