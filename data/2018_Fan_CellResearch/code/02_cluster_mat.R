source('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/resource/01_function.R')
mat = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Fan_CellResearch/proc/expr.rds')
mat = log2(mat+1)
saveRDS(mat,'/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Fan_CellResearch/proc/expr.rds')
getMatrixCluster(data.path = '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Fan_CellResearch/proc/expr.rds', result.path = '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Fan_CellResearch/proc/', platform = 'none')

