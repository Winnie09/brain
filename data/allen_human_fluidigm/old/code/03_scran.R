suppressMessages(library(scran))
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/allenBrain/data/')
cnt = readRDS('./proc/matrix/count.rds')
sce <- SingleCellExperiment(list(counts=cnt))
sce <- computeSumFactors(sce,BPPARAM=MulticoreParam(workers = 10))
sf <- sizeFactors(sce)
saveRDS(sf,'./proc/scran/sizefactor.rds')
normmat <- sweep(cnt,2,sf,'/')
# normmat <- normmat[rowMeans(normmat > 0) > 0.01,]
dir.create('./proc/scran/',showWarnings = F, recursive = T)
saveRDS(normmat,file='./proc/matrix/scran.rds')  

