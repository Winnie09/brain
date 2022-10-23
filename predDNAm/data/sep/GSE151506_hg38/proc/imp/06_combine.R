rm(list=ls())
ddir <- '/home/whou10/scratch16/whou10/brain/predDNAm/data/sep/GSE151506_hg38/proc/imp/indimp/scRNA/'
rdir <- '/home/whou10/scratch16/whou10/brain/predDNAm/data/sep/GSE151506_hg38/proc/imp/'
af <- list.files(ddir)
af <- list.files(ddir)
ad <- lapply(af, function(f){
  d <- readRDS(paste0(ddir, f))
})
ad <- do.call(cbind, ad)
saveRDS(ad, paste0(rdir, 'scRNA_paired_saver_tpm.rds'))
saveRDS(log2(ad + 1), paste0(rdir, 'scRNA_paired_saver_log2tpm.rds'))
