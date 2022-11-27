m <- readRDS('/home/whou10/scratch16/whou10/brain/predDNAm/data/sep/GSE151506_hg38/proc/scRRBS/promoter_up1k_down1k.rds')
rdir <- '/home/whou10/scratch16/whou10/brain/predDNAm/knn/umap_dist_vs_promoter_DNAme_dist/res_dlist_single/'
dir.create(rdir, showWarnings = F, recursive = T)
library(parallel)
cut = as.numeric(commandArgs(trailingOnly = T)[[1]][1])
print(cut)
res <- lapply((1+(cut-1)*10): min((cut*10), ncol(m)), function(i){
  print(i)
  tmpv <- rep(NA, ncol(m))
  a <- mclapply((i+1):(ncol(m)), function(j){
    mean(m[,i] == m[,j], na.rm = T)
  }, mc.cores = 24)
  tmpv[(i+1):ncol(m)] <- unlist(a)
  names(tmpv) <- colnames(m)
  tmpv
})
names(res) <- colnames(m)[(1+(cut-1)*10): min((cut*10), ncol(m))]
saveRDS(res, paste0(rdir, 'dist_DNAme_proportion_agreement_list_',
                    (1+(cut-1)*10), '_', 
                    min((cut*10), ncol(m)), '.rds'))

