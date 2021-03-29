rm(list = ls())
library(here)
setwd(here())
source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R')
ddir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/mouseGBM/pseudotime/res/testvar_saver_pm/'

alltraj = list.files(ddir)
alltraj = alltraj[grepl('_', alltraj)]

for (traj in alltraj){
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
      Res = readRDS(paste0(rdir, 'res.rds'))
      stat = Res$statistics
      stat <- stat[order(stat[,1], -stat[,3]), ]
      diffgene <- rownames(stat[stat[, grep('^fdr.*overall$', colnames(stat))] < 0.05,])
      
      
      str(diffgene)
      if (length(diffgene) > 0){
        write.csv(stat[diffgene, , drop=F], paste0(pdir, 'differential_genes.csv'))  
        
        Res$populationFit <- getPopulationFit(Res, gene = diffgene, type = 'variable')
        
        ## -----------
        ## clustering
        ## -----------
        
        Res$covariateGroupDiff <- getCovariateGroupDiff(testobj = Res, gene = diffgene, reverse = F)
        
        DEGType <- getDEGType(Res)
        sink(paste0(pdir, 'DEGType_table.txt'))
        table(DEGType)
        sink()
        
        if (length(diffgene) > 20){
          clu <- clusterGene(Res, gene = names(DEGType)[!DEGType %in% c('nonDEG', 'meanSig')], type = 'variable', k=3)
          table(clu)
          # clu2 <- rep(6, sum(DEGType == 'meanSig'))
          # if (length(clu2) > 0) {
          #   names(clu2) <- names(DEGType)[DEGType %in% c('meanSig')]
          #   clu <- c(clu, clu2)
          # }
          design = Res$design
          cellanno = Res$cellanno
          meandiff <- sapply(c(0,1), function(i){
            as <- rownames(design[design[,2]==i, ])
            rowMeans(Res$expr[names(DEGType)[DEGType == 'meanSig'], cellanno[cellanno[,2] %in% as,1]])
          })
          
          large0 <- rownames(meandiff)[meandiff[,1] >= meandiff[,2]]
          large1 <- rownames(meandiff)[meandiff[,1] < meandiff[,2]]
          
          clu2 <- rep(6, length(large0))
          names(clu2) <- large0
          clu3 <- rep(7, length(large1))
          names(clu3) <- large1
          clu = c(clu, clu2, clu3)
          
        } else {
          clu <- rep(1, length(diffgene))
          names(clu) <- diffgene
        }
        
        
        Res$cluster <- clu
        saveRDS(Res, paste0(pdir, 'numeric_res_with_clu.rds'))
        ## --------------
        ## save diff gene
        ## --------------
        
        res <- data.frame(gene = diffgene, stat[diffgene, ,drop=F], cluster = Res$cluster[diffgene], stringsAsFactors = F)
        res <- cbind(res, DEGType = DEGType[rownames(res)])
        write.csv(res, paste0(pdir, 'differential_genes_full.csv'))
        
        ## -----------
        ## GO analysis
        ## -----------
        if (length(diffgene) >= 20){
          goRes <- GOEnrich(testobj = Res, type = 'variable')
          saveRDS(goRes, paste0(pdir, '/goRes.rds'))
          
          nn <- sapply(1:length(goRes), function(i){
            tmp <- goRes[[i]]
            tmp <- tmp[tmp[, 'FDR'] < 0.05, , drop=F]
            print(str(tmp))
            if (nrow(tmp) > 0){
              write.csv(tmp, paste0(pdir, 'cluster', i, '_GO.csv'))  
              return(nrow(tmp))
            } else {
              return(0)  
            }
          
          })
          if (sum(nn) > 0){
            pdf(paste0(pdir, '/hm_GO_term.pdf'), width = 7.2, height = 3.5)
            print(plotGOEnrich(goRes))
            dev.off()
          }
            
        }
          
        

        
        ## -----------------------
        ## plotClusterMeanAndDiff
        ## -----------------------
        pdf(paste0(pdir, 'cluster_mean_and_diff.pdf'), width = 3.8, height = 7.5)
        print(plotClusterMeanAndDiff(Res, cluster = Res$cluster))
        dev.off()
        
        # --------------------------------------
        # compare original and fitted expression
        # --------------------------------------
        if (length(diffgene) < 100){
          png(paste0(pdir, 'DiffFitHm_rownames.png'),width = 12000,height = 10000,res = 300)
          plotDiffFitHm(Res, type='variable', subsampleCell = T,  numSubsampleCell = min(ncol(Res$covariateGroupDiff), 5e3), showRowName = T, cellWidthTotal = 300, cellHeightTotal = length(Res$cluster) * 10)
          dev.off()
        } else {
          png(paste0(pdir, 'DiffFitHm.png'),width = 4000,height = 2200,res = 200)
          plotDiffFitHm(Res, type = 'variable', cellWidthTotal = 200, cellHeightTotal = 300)
          dev.off()  
        }
          
              
        ## ----------
        ## plot DEG
        ## ----------
        DEGType <- DEGType[diffgene]
        id <- sort(sample(1:ncol(Res$populationFit[[1]]), ncol(Res$expr)))
        Res$populationFit[[1]] <- Res$populationFit[[1]][, id]
        Res$populationFit[[2]] <- Res$populationFit[[2]][, id]
        
        for (i in unique(DEGType)){  ## debug -- ok!!
          print(i)
          gene <- names(DEGType)[DEGType == i]
          png(paste0(pdir, 'diffgene_sampleFit_', i, '.png'), width = 4000, height = 2500, res = 200)
          print(plotGene(Res, gene = gene[1:min(length(gene), 25)], plot.point = T, point.size = 0.1, variable = 'location'))
          dev.off()
        }
        
        for (i in unique(DEGType)){
          print(i)
          gene <- names(DEGType)[DEGType == i]
          png(paste0(pdir, 'diffgene_groupFit_', i, '.png'), width = 2500, height = 2500, res = 200)
          print(plotGenePopulation(testobj = Res, type = 'variable', gene = gene[1:min(length(gene), 100)]))
          dev.off()
        }
        
        for (i in 1:max(Res$cluster)){
          print(i)
          gene <- rownames(res)[res$cluster == i]
          png(paste0(pdir, 'diffgene_groupFit_cluster', i, '.png'), width = 2500, height = 2500, res = 200)
          print(plotGenePopulation(testobj = Res, type = 'variable', gene = gene[1:min(length(gene), 100)]))
          dev.off()
        } 
      }
    }
  }
}






