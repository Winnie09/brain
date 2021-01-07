library(here)
setwd(here())
gbm.meta <- readRDS('atlasGBM/GBMonly/data/36nonNormal_combined_meta.rds')
atlas.meta <- readRDS('atlasGBM/humanatlas/data/meta_allcell.rds')
colnames(gbm.meta)[colnames(gbm.meta) == 'Age'] <- 'age'
colnames(gbm.meta)[colnames(gbm.meta) == 'Sex'] <- 'sex'
colnames(gbm.meta)[colnames(gbm.meta) == 'Grade'] <- 'grade'
colnames(gbm.meta)[colnames(gbm.meta) == 'Location'] <- 'location'
colnames(gbm.meta)[colnames(gbm.meta) == 'Treatment'] <- 'treatment'
colnames(atlas.meta)[colnames(atlas.meta) == 'gender'] <- 'sex'
colnames(atlas.meta)[colnames(atlas.meta) == 'tumor.grade'] <- 'Tumor.Grade'
colnames(atlas.meta)[colnames(atlas.meta) == 'Treatment'] <- 'treatment'
atlas.meta <- atlas.meta[, -which(colnames(atlas.meta) == 'sample.id')] 
colnames(atlas.meta)[colnames(atlas.meta) == 'donor'] <- 'Sample.ID'

meta <- reshape::merge_recurse(list(atlas.meta, gbm.meta))
saveRDS(meta, 'atlasGBM/humanGBM/data/36nonNormalGBM_humanatlas_combined_meta.rds')

gbm <- readRDS('atlasGBM/GBMonly/data/36nonNormal_combined_log2norm_dlist.rds')
atlas <- readRDS('atlasGBM/humanatlas/data/dlist.rds')
dlist <- c(atlas, gbm)
intgn <- sapply(dlist,row.names)
intgn <- table(unlist(intgn))
intgn <- names(intgn)[intgn==length(dlist)]
for (i in names(dlist)){
  dlist[[i]] <- dlist[[i]][intgn, ]
}
saveRDS(dlist, 'atlasGBM/humanGBM/data/36nonNormalGBM_humanatlas_combined_log2norm_dlist.rds')

mat <- do.call(cbind, dlist)
saveRDS(mat, 'atlasGBM/humanGBM/data/36nonNormalGBM_humanatlas_combined_log2norm_mat.rds')

