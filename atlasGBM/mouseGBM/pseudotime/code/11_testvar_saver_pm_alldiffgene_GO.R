rm(list = ls())
library(here)
setwd(here())
source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R')
ddir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/mouseGBM/pseudotime/res/testvar_saver_pm/'

alltraj = list.files(ddir)
alltraj = alltraj[grepl('_', alltraj)]

for (traj in alltraj){
  # traj = as.character(commandArgs(trailingOnly = T)[1])
  # comparison = as.character(commandArgs(trailingOnly = T)[2])
  print(traj)
  allcomparison = list.files(paste0(ddir, traj))
  allcomparison = allcomparison[grepl('_', allcomparison)]
  for (comparison in allcomparison){
    print(comparison)
    rdir <- paste0(ddir, traj, '/', comparison, '/')
    
    if (file.exists(paste0(rdir, '/res.rds'))){
      pdir <- paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/mouseGBM/pseudotime/plot/testvar_saver_pm/', traj, '/', comparison, '/')
      dir.create(pdir, showWarnings = F, recursive = T)
      # ---------------------
      # prepare data and test
      # ---------------------
      if (file.exists(paste0(rdir, 'numeric_res_with_clu.rds'))){
        Res <-  readRDS(paste0(rdir, 'numeric_res_with_clu.rds'))
      } else if (file.exists(paste0(pdir, 'numeric_res_with_clu.rds'))){
        Res <- readRDS(paste0(pdir, 'numeric_res_with_clu.rds'))
      }
      
      stat = Res$statistics
      stat <- stat[order(stat[,1], -stat[,3]), ]
      diffgene <- rownames(stat[stat[, grep('^fdr.*overall$', colnames(stat))] < 0.05,])
      
      ## -----------
      ## GO analysis
      ## -----------
      
      if (length(diffgene) >= 5){
        goRes <- GOEnrich(testobj = Res, type = 'variable', use.clusters = F)
        saveRDS(goRes, paste0(pdir, 'goRes_alldiffgene.rds'))
        
        tmp <- goRes[[1]]
        tmp <- tmp[tmp[, 'FDR'] < 0.05, , drop=F]
        print(str(tmp))
        if (nrow(tmp) > 0){
          write.csv(tmp, paste0(pdir, 'allgene_GO.csv'))  
          pdf(paste0(pdir, 'hm_GO_term_alldiffgene.pdf'), width = 7.2, height = 3.5)
          print(plotGOEnrich(goRes))
          dev.off() 
        }
      }
    }
  }
}
