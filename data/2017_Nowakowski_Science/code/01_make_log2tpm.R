library(data.table)
d = fread('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2017_Nowakowski_Science/raw/exprMatrix.tsv.gz', data.table = F)
gn = d[-1,1]
d = as.matrix(d[,-1])
d = d[-1,]
# d = matrix(as.numeric(d),nrow=length(gn))
rownames(d) = gn
colnames(d) = paste0('2017_Nowakowski_Science:',colnames(d))
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/function.R')
a = log2TPM_from_fludigm_count(d)
saveRDS(a, '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2017_Nowakowski_Science/proc/log2tpm.rds')

