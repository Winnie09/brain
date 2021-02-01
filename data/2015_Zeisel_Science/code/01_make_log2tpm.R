library(data.table)
d = fread('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2015_Zeisel_Science/raw/GSE60361_C1-3005-Expression.txt', data.table = F)
gn = d[-1,1]
d = as.matrix(d[,-1])
d = d[-1,]
rownames(d) = gn
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/myfunc/TPM_from_fludigm_count.R')
a = TPM_from_fludigm_count(d, species = 'mouse', referenceGenome = 'grcm38', lib.factor = 1e6)
range(a)
colnames(a) = paste0('2015_Zeisel_Science:', colnames(a))
saveRDS(a, '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2015_Zeisel_Science/proc/tpm.rds')

source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/myfunc/log2TPM_from_fludigm_count.R')
a = log2TPM_from_fludigm_count(d, species = 'mouse', referenceGenome = 'grcm38', lib.factor = 1e6)
colnames(a) = paste0('2015_Zeisel_Science:', colnames(a))
saveRDS(a, '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2015_Zeisel_Science/proc/log2tpm.rds')


