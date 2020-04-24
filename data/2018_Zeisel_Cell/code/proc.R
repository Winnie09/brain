mat <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Zeisel_Cell/proc/mat.rds')
meta <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Zeisel_Cell/raw/meta/cell.rds')
cn <- paste0('2018_Zeisel:cell',1:ncol(mat))
meta[,2] <- as.character(meta[,2])
meta[,3] <- as.character(meta[,3])
meta[,4] <- as.character(meta[,4])
meta[,1] <- cn
colnames(mat) <- cn
meta <- data.frame(cell=cn,Location=meta[,2],celltype=meta[,3],time=meta[,4],stringsAsFactors = F)
saveRDS(mat,'/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Zeisel_Cell/proc/mat.rds')
saveRDS(meta,'/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Zeisel_Cell/proc/meta.rds')
