num <- as.numeric(commandArgs(trailingOnly = TRUE)[1])
print(num)
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/myfunc/01_function.R')
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/')
dlist <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/humanatlas/data/dlist_add_2020_Guo_Nature.rds')
i = names(dlist)[num]
print(i)
mat = dlist[[i]]
rdir <- pdir <- paste0(i, '/pca/')
dir.create(rdir, showWarnings = FALSE, recursive = TRUE)
pr <- PCA(genebycellmat = mat, save.pca = TRUE, plot.statistics=FALSE, plot.dir = pdir, result.dir = rdir, PC_for_columns = TRUE, findVariableGenes = TRUE, findVariableGenesMethod = 'lm', maxVariableGenes = NULL, numPC = 10, smoothFittingMethod = 'loess')

