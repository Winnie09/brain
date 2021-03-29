rm(list = ls())
library(here)
setwd(here())
source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R')
alltraj = list.files('atlasGBM/mouseGBM/pseudotime/res/testvar_saver_pm/')
alltraj = alltraj[grepl('_', alltraj)]
res <- sapply(alltraj, function(traj){
  print(traj)
  allcomparison = list.files(paste0('atlasGBM/mouseGBM/pseudotime/res/testvar_saver_pm/', traj))
  allcomparison <- allcomparison[grepl('_', allcomparison)]
  
  sapply(allcomparison, function(comparison){
    print(comparison)
    rdir <- paste0('atlasGBM/mouseGBM/pseudotime/res/testvar_saver_pm/', traj, '/', comparison)
    
    if (file.exists(paste0(rdir, '/res.rds'))){
      pdir <- paste0('atlasGBM/mouseGBM/pseudotime/plot/testvar_saver_pm/', traj, '/', comparison)
      dir.create(pdir, showWarnings = F, recursive = T)
      # ---------------------
      # prepare data and test
      # ---------------------
      Res = readRDS(paste0(rdir, '/res.rds'))
      stat = Res$statistics
      stat <- stat[order(stat[,1]), ]
      diffgene <- rownames(stat[stat[, grep('^fdr.*overall$', colnames(stat))] < 0.05,])
      return(length(diffgene))
    }else{
      return(NA)
    }
  })  
})

res <- do.call(cbind, res)
dir.create('atlasGBM/mouseGBM/pseudotime/res/testvar_saver_pm/perf/', showWarnings = F, recursive = T)
saveRDS(res, 'atlasGBM/mouseGBM/pseudotime/res/testvar_saver_pm/perf/number_of_diffgene.rds')
write.csv(res, 'atlasGBM/mouseGBM/pseudotime/res/testvar_saver_pm/perf/number_of_diffgene.csv')
