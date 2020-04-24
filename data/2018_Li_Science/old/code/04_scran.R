suppressMessages(library(scran))
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Li_Science/')
af = list.files('./processed/qc/')
for (f in af){
  cnt = readRDS(paste0('./processed/qc/',f))
  sce <- SingleCellExperiment(list(counts=cnt))
  sce <- computeSumFactors(sce,BPPARAM=MulticoreParam(workers = 10))
  sf <- sizeFactors(sce)
  normmat <- sweep(cnt,2,sf,'/')
  # normmat <- normmat[rowMeans(normmat > 0) > 0.01,]
  dir.create('./processed/scran/',showWarnings = F, recursive = T)
  saveRDS(normmat,file=paste0('./processed/scran/',f))  
}

