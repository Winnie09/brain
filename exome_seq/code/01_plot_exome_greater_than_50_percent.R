af <- list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/exome_seq/data/', pattern = '.csv')
tb <- read.csv(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/exome_seq/data/', af[1]), as.is = T)

for (f in af[2:length(af)]){
  print(f)
  tmp <- read.csv(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/exome_seq/data/', f), as.is = T)
  tb <- rbind(tb, tmp)
}

tb <- tb[grepl('GBM', tb[,2]), ]

tb[,2] <- sapply(tb[,2], function(i) substr(i, 1, 6))

mat <- matrix(0, nrow = length(unique(tb[,2])), ncol = 44)
rownames(mat) <- unique(tb[,2])
colnames(mat) <- c(paste0(paste0('chr', 1:22), '+'), paste0(paste0('chr', 1:22), '-'))

for (p in unique(tb[,2])){
  print(p)
  tmp = tb[tb[,2]==p, ]
  for (i in 1:22){
    print(i)
    v <- which(tmp[, paste0('chr', i)] > 50)
    if (length(v) > 0) {
      if (sum(grepl('Amplification', tmp[v,1]) & grepl('3', tmp[v,1])) >= 1) mat[p, paste0('chr', i, '+')] <- 3
      if (sum(grepl('Amplification', tmp[v,1]) & grepl('4', tmp[v,1])) >= 1) mat[p, paste0('chr', i, '+')] <- 4
      if (sum(grepl('Amplification', tmp[v,1]) & grepl('5', tmp[v,1])) >= 1) mat[p, paste0('chr', i, '+')] <- 5
      if (sum(grepl('Amplification', tmp[v,1]) & grepl('6', tmp[v,1])) >= 1) mat[p, paste0('chr', i, '+')] <- 6
      if (sum(grepl('Amplification', tmp[v,1]) & grepl('7', tmp[v,1])) >= 1) mat[p, paste0('chr', i, '+')] <- 7
      if (grepl('Deletion: Homozygous', tmp[v,1])) mat[p, paste0('chr', i, '-')] <- -3
      if (grepl('Deletion: Loss of Heterozygosity', tmp[v,1])) mat[p, paste0('chr', i, '-')] <- -6
    }  
  }  
}

library(pheatmap)
pdir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/exome_seq/plot/'
pdf(paste0(pdir, 'exome_greater_than_50_percent.pdf'), height = 4, width = 7.2)
pheatmap(mat, cluster_rows = F, cluster_cols = F)
dev.off()

saveRDS(mat, paste0(pdir, 'exome_greater_than_50_percent.rds'))

