rm(list=ls())
library(here)
setwd(here())
source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R')
ddir = '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/trajectory/res/twogroup/'
dgList <- sapply(1:6, function(i){
  print(i)
  comp = sub('.rds', '', list.files(ddir, pattern = '.rds'))[i]
  ## "Frontal_Occipital": frontal 1, occipital 0 (higher expression)
  Res <- readRDS(paste0(ddir, comp, '.rds'))
  print(names(Res))
  str(Res$llr.overall)
  stat = Res$statistics
  diffgene <- rownames(stat[stat[, grep('^fdr.*overall$', colnames(stat))] < 0.05,])
})
saveRDS(dgList, '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/trajectory/res/twogroup/perf/dgList.rds')
