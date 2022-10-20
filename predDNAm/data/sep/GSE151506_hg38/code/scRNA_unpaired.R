rm(list=ls())
ddir <- '/home/whou10/scratch16/whou10/brain/predDNAm/data/sep/GSE151506_hg38/raw/scRNA/'
rdir <- '/home/whou10/scratch16/whou10/brain/predDNAm/data/sep/GSE151506_hg38/proc/scRNA/'
af = list.files(ddir)
alld <- lapply(af, function(f){
  print(f)
  d <- read.table(paste0(ddir, f), sep = '\t')
  colnames(d) <- d[1,]
  d <- d[-1,]
  rownames(d) <- d[,1]
  d <- d[,-1]
  dn <- dimnames(d)
  d2 <- matrix(as.numeric(as.matrix(d)), nrow = nrow(d))
  dimnames(d2) <- dn
  d2
})

alld2 <- do.call(cbind, alld)
saveRDS(alld2, paste0(rdir, 'tpm_unpaired.rds'))
