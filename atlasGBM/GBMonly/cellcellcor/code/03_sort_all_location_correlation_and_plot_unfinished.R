library(RColorBrewer)
library(pheatmap)
library(gplots)
library(reshape2)
library(ggplot2)
library(gplots)
library(parallel)
order.loc <- c('Frontal', 'Occipital', 'Parietal', 'Temporal')
co1 <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cellcellcor/res/cellorder_Frontal.rds')) ##### add
co2 <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cellcellcor/res/cellorder_Occipital.rds')) ##### add
co3 <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cellcellcor/res/cellorder_Parietal.rds')) ##### add
co4 <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cellcellcor/res/cellorder_Temporal.rds')) ##### add
co <- c(co1, co2, co3, co4)

## ============
mat <- lapply(order.loc, function(loc1){
  # loc1 = 'Frontal'
  print(loc1)
  cellcor <- mclapply(setdiff(order.loc, loc1), function(loc2){
    print(loc2)
    fn <- paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cellcellcor/res/correlationmatrix_', loc1, '_', loc2, '_centered.rds')
    if (file.exists(fn)){
      tmp = readRDS(fn)
    } else {
      tmp = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cellcellcor/res/correlationmatrix_', loc2, '_', loc1, '_centered.rds'))
      tmp <- t(tmp)
    }
  }, mc.cores = 3)
  names(cellcor) <- setdiff(order.loc, loc1)
  cellcor[[loc1]] <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cellcellcor/res/correlationmatrix_', loc1, '_centered.rds'))
  cellcor <- cellcor[order.loc]
  cellcor <- do.call(cbind, cellcor)
  print(str(cellcor))  
  cellcor <- cellcor[, co]
  # cellcor <- cellcor[1:10, 1:100] ####  
})
saveRDS(mat, '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cellcellcor/res/correlationmatrix_all_locations_list.rds')

mat <- do.call(rbind, mat)
mat <- mat[co, ]
identical(colnames(mat), co)
allp <- sub('_.*', '', co)

png(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cellcellcor/plot/cellcor_all_locations_centered.png'), width = 8000, height = 8000, res = 100)
heatmap(mat, scale = 'none', Rowv = NA, Colv = NA, labRow = NA, labCol = NA,col=colorRampPalette(c('blue3','blue1','skyblue','cyan','green','lightgreen','yellow','orange','pink','red1','red3'))(100),add.expr=eval({abline(v=cumsum(rle(allp)$lengths), lwd=5);abline(h=cumsum(rle(allp)$lengths), lwd=5)}))
dev.off()


