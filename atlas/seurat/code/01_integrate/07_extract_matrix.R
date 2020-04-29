key = as.character(commandArgs(trailingOnly = T)[[1]])
setwd(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlas/result/',key))
d = readRDS('./integrated.rds')
mat1 = d@assays$integrated@data
mat2 = d@assays$RNA@data
saveRDS(mat1, './mat_after_Seurat.rds')
saveRDS(mat2, './mat_before_Seurat.rds')
