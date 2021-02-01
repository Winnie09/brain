library(data.table)
d = fread('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2016_Tasic_NatNeuro/raw/GSE71585_RefSeq_TPM.csv.gz', data.table = F)
gn = d[-1,1]
d = as.matrix(d[,-1])
d = d[-1,]
rownames(d) = gn
str(d)
range(d)
saveRDS(d,'/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2016_Tasic_NatNeuro/proc/tpm.rds')
d = log2(d + 1)
colnames(d) = paste0('2016_Tasic_NatNeuro:',colnames(d))
str(d)
range(d)
saveRDS(d,'/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2016_Tasic_NatNeuro/proc/log2tpm.rds')

