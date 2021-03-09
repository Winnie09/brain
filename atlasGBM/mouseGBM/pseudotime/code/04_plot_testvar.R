source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R')
rdir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/mouseGBM/pseudotime/res/testvar/'
pdir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/mouseGBM/pseudotime/plot/testvar/'

af <- list.files(rdir)
af <- af[!grepl('pt', af)]
af <- af[grepl('_', af)]
f = af[1]

for (f in rev(af)){
  print(f)
  dir.create(paste0(pdir, f), recursive = T, showWarnings = F)
  if (file.exists(paste0(rdir, f, '/res.rds'))){
    Res = readRDS(paste0(rdir, f, '/res.rds'))
    statistics = Res$statistics
    deg <- rownames(statistics)[statistics$both.fdr < 0.05]
    res <- data.frame(gene = deg, DEGType = getDEGType(Res)[deg], stringsAsFactors = F)
    res <- cbind(res, statistics[deg, ])
    res <- res[order(res$both.fdr, -res$both.fc), ]
    deg <- rownames(res)
    write.csv(res, paste0(rdir, f, '/DEG.csv'))
    saveRDS(res, paste0(rdir, f, '/DEG.rds'))
    
    if (length(deg) == 1){
        png(paste0(pdir, f, '/diffgene_groupFit_', deg, '.png'), width = 800, height = 800, res = 200)
        print(plotGenePopulation(testobj = Res, type = 'variable', gene = deg))
        dev.off()
        
        Res$populationFit <- getPopulationFit(Res, gene = deg, type = 'variable')
        png(paste0(pdir, f, '/diffgene_groupDiff_', deg, '.png'), width = 500, height = 500, res = 200)
        print(plotClusterDiff(testobj = Res, gene = deg, each = TRUE, sep = ':.*'))
        dev.off()
        
    }
    
    if (length(deg) > 1){
      ## --------------
      ## population fit
      ## --------------
      Res$populationFit <- getPopulationFit(Res, gene = deg, type = 'variable')
      
      ## -----------
      ## clustering
      ## -----------
      Res$covariateGroupDiff <- getCovariateGroupDiff(testobj = Res, gene = rownames(statistics[statistics$both.fdr < 0.05,]))
      
      if (length(deg) > 10){
        Res$cluster <- clusterGene(Res, gene = deg, type = 'variable', k=3)
      } else {
        Res$cluster <- clusterGene(Res, gene = deg, type = 'variable', k=1)
      }
        
      
      ## ----------------
      ## plotClusterDiff
      ## ----------------
      pdf(paste0(pdir, f, '/cluster_diff.pdf'), width = 3, height = 2)
      plotClusterDiff(testobj=Res, gene = deg)
      dev.off()

      ## ---------------
      ## plotClusterMean
      ## ----------------
      pdf(paste0(pdir, f, '/cluster_mean.pdf'), width = 5.5, height = 2.2)
      plotClusterMean(testobj=Res, cluster = Res$cluster, type = 'variable')
      dev.off()

      ## -----------
      ## GO analysis
      ## -----------
      goRes <- GOEnrich(testobj = Res, type = 'variable', version = 3)
      nn <- sapply(1:length(goRes), function(i){
        tmp <- goRes[[i]]
        tmp <- tmp[tmp[, 'FDR'] < 0.05, ]
        write.csv(tmp, paste0(rdir, f, '/cluster', i, '_GO.csv'))
        print(str(tmp))
        return(0)
      })

      ## -----------------------
      ## plotClusterMeanAndDiff
      ## -----------------------
      pdf(paste0(pdir, f, '/cluster_mean_and_diff.pdf'), width = 4, height = 5.5)
      print(plotClusterMeanAndDiff(Res, cluster = Res$cluster))
      dev.off()

      res$cluster <- Res$cluster[rownames(res)]
      for (i in 1:max(Res$cluster)){
        print(i)
        gene <- rownames(res[res$cluster == i, ])
        png(paste0(pdir, f, '/diffgene_groupFit_cluster', i, '.png'), width = 1500, height = 1500, res = 200)
        print(plotGenePopulation(testobj = Res, type = 'variable', gene = gene[1:min(length(gene), 100)]))
        dev.off()
      }

      for (i in 1:max(Res$cluster)){
        print(i)
        gene <- rownames(res[res$cluster == i, ])
        png(paste0(pdir, f, '/diffgene_groupDiff_cluster', i, '.png'), width = 1500, height = 1500, res = 200)
        print(plotClusterDiff(testobj = Res, gene = gene[1:min(length(gene), 100)], each = TRUE, sep = ':.*'))
        dev.off()
      }

    
      # --------------------------------------
      # compare original and fitted expression
      # --------------------------------------
      png(paste0(pdir, f, '/fitHm.png'),width = 4000,height = 2500,res = 300)
      plotFitHm(Res, type = 'variable',showRowName = T,subsampleCell=F)
      dev.off()  
    }
  }
}
  
