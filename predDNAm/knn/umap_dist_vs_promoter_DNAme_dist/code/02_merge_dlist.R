rm(list=ls())
ddir <- '/home/whou10/scratch16/whou10/brain/predDNAm/knn/umap_dist_vs_promoter_DNAme_dist/res_dlist_single/'
af <- list.files(ddir, pattern = 'dist_DNAme_proportion_agreement_list_.*_.*rds')
cut <- sapply(af, function(i){
  as.numeric(sub('_.*','', sub('dist_DNAme_proportion_agreement_list_', '', i)))
})
names(cut) <- af
cut <- sort(cut)
dlist <- lapply(1:length(cut), function(i){
  d <- readRDS(paste0(ddir, names(cut)[i]))
  d <- do.call(rbind, d)
  if (ncol(d) > 1313) d <- d[, 1:1313]
  return(d)
})
d3 <- do.call(rbind, dlist)
dm <- dimnames(d3)
d4 <- matrix(as.numeric(d3), nrow = nrow(d3))
dimnames(d4) <- dm
saveRDS(d4, '/home/whou10/scratch16/whou10/brain/predDNAm/knn/umap_dist_vs_promoter_DNAme_dist/res/DNAme_proportion_agreement.rds')
