library(destiny)
# study <- '2018_Li_Science'
study <- as.character(commandArgs(trailingOnly = TRUE)[[1]])
print(study)
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/myfunc/01_function.R')
ddir <- paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/', study, '/pca/')
pdir <- rdir <- paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/', study)
pr <- readRDS(paste0(ddir, 'pr.rds'))
UMAP(samplebyfeature_mat = pr, save.umap = TRUE, result.dir = rdir)
