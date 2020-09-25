library(data.table)
d <- fread('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2019_loo_natcom/raw/GSE123335_E14_combined_matrix.txt.gz',data.table=F)
row.names(d) <- d[,1]
d <- as.matrix(d[,-1])

dm <- fread('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2019_loo_natcom/raw/GSE123335_E14_combined_matrix_ClusterAnnotations.txt.gz',data.table=F)
dm[,1] <- gsub('\\.','-',dm[,1])
dm[,2] <- sub(' .*','',dm[,2])

m <- data.frame(Cell=colnames(d),Celltype=dm[match(colnames(d),dm[,1]),2])
m$Location <- 'cerebral cortex'
m$Time <- 'E14.5'
saveRDS(d,file='/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2019_loo_natcom/proc/expr_E14.rds')
saveRDS(m,file='/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2019_loo_natcom/proc/meta_E14.rds')

d <- fread('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2019_loo_natcom/raw/GSE123335_P0_combined_matrix.txt.gz',data.table=F)
row.names(d) <- d[,1]
d <- as.matrix(d[,-1])

dm <- fread('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2019_loo_natcom/raw/GSE123335_P0_combined_matrix_ClusterAnnotations.txt.gz',data.table=F)
dm[,1] <- gsub('\\.','-',dm[,1])
dm[,2] <- sub(' .*','',dm[,2])

m <- data.frame(Cell=colnames(d),Celltype=dm[match(colnames(d),dm[,1]),2])
m$Location <- 'cerebral cortex'
m$Time <- 'P0'
saveRDS(d,file='/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2019_loo_natcom/proc/expr_P0.rds')
saveRDS(m,file='/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2019_loo_natcom/proc/meta_P0.rds')

