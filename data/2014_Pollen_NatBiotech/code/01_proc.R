a = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2014_Pollen_NatBiotech/raw/pollen.rds')
m=assays(a)$normcounts
saveRDS(m, '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2014_Pollen_NatBiotech/proc/tpm.rds')

log2m = log2(m + 1)
saveRDS(log2m, '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2014_Pollen_NatBiotech/proc/log2tpm.rds')

ct = sapply(colnames(m), function(i) strsplit(i, '_')[[1]][2])
meta = data.frame(cell = names(ct), celltype = ct, stringsAsFactors = FALSE)
saveRDS(meta, '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2014_Pollen_NatBiotech/proc/meta.rds')
