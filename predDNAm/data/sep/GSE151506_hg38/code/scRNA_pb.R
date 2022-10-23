## non-imputed values
d2 <- readRDS('/home/whou10/scratch16/whou10/brain/predDNAm/data/sep/GSE151506_hg38/proc/scRNA/tpm_paired.rds')
d <- log2(d + 1)
s <- sub(':.*','',colnames(d))
m <- rowsum(t(d),s)
m <- m/as.vector(table(s)[rownames(m)])
m <- t(m)
saveRDS(m,file='/home/whou10/scratch16/whou10/brain/predDNAm/data/sep/GSE151506_hg38/proc/scRNA/log2tpm_pb.rds')

## imputed values 
d <- readRDS('/home/whou10/scratch16/whou10/brain/predDNAm/data/sep/GSE151506_hg38/proc/imp/scRNA_paired_saver_log2tpm.rds')
s <- sub(':.*','',colnames(d))
m <- rowsum(t(d),s)
m <- m/as.vector(table(s)[rownames(m)])
m <- t(m)
saveRDS(m,file='/home/whou10/scratch16/whou10/brain/predDNAm/data/sep/GSE151506_hg38/proc/scRNA/saver_log2tpm_pb.rds')

