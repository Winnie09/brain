setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/')     
rdir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/humanatlas/data/'
mat <- readRDS(paste0(rdir, 'combine_mat_add_2020_Guo_Nature.rds'))
meta.new2 <- readRDS(paste0(rdir, 'meta_allcell.rds'))

#### add HCL meta
meta <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2020_Guo_Nature/meta/meta_sub.rds')
meta$cell <- paste0('2020_Guo_Nature;', meta$cell)
meta.na <- matrix(NA, nrow = nrow(meta), ncol = (ncol(meta.new2) - ncol(meta)))
dimnames(meta.na) <- list(rownames(meta), setdiff(colnames(meta.new2), colnames(meta)))
meta <- cbind(meta, meta.na)
meta <- meta[, colnames(meta.new2)]
rownames(meta) <- meta$cell

c.s <- setdiff(colnames(mat), rownames(meta.new2))
meta <- meta[c.s, ]

meta.new3 <- rbind(meta.new2, meta)
identical(meta.new3$cell, colnames(mat))
saveRDS(meta.new3, paste0(rdir, 'meta_allcell_add_2020_Guo_Nature.rds'))

