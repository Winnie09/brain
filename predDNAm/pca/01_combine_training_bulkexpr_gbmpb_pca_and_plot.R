## read in pseudobulk
m <- readRDS(file = '/home/whou10/scratch16/whou10/brain/atlasGBM/GBMonly/pseudobulk/data/exprpb_sample_level.rds')
train.pb <- readRDS('/home/whou10/scratch16/whou10/brain/predDNAm/data/combine/res/combinege.rds')

## combinei on intersect genes
int = intersect(rownames(m), rownames(train.pb))
gbm.pb = m[int, ]
colnames(gbm.pb) <- paste0('JHU:', colnames(gbm.pb))
train.pb = train.pb[int, ]
combine <- cbind(gbm.pb, train.pb)
stu <- sub(':.*', '', colnames(combine))
saveRDS(combine, '/home/whou10/scratch16/whou10/brain/predDNAm/pca/combine_training_gbm_pseudobulk.rds')

combine.bak = combine

## retain highly expressed genes
library(preprocessCore)
dn <- dimnames(combine)
combine <- normalize.quantiles(combine)
dimnames(combine) <- dn
combine <- combine[rowMeans(combine > 1) > 0.1,  ]

## pca
source('/home/whou10/scratch16/whou10/resource/myfunc/01_function.R')
PCA(genebycellmat = combine, save.pca = T, plot.statistics=F, plot.dir = '/home/whou10/scratch16/whou10/brain/predDNAm/pca/', result.dir = '/home/whou10/scratch16/whou10/brain/predDNAm/pca/', PC_for_columns = TRUE, findVariableGenes = TRUE, findVariableGenesMethod = 'lm', maxVariableGenes = 2e3, numPC = 10, smoothFittingMethod = 'loess')

pr  <- readRDS('/home/whou10/scratch16/whou10/brain/predDNAm/pca/pr.rds')
stu <- sub(':.*', '', rownames(pr))

## plot
pd = data.frame(PC1 = pr[,1], PC2 = pr[,2], stringsAsFactors = F, study = stu)
rownames(pd) <- rownames(pr)
library(ggplot2)
pdf('/home/whou10/scratch16/whou10/brain/predDNAm/pca/pca_study.pdf', height = 4, width = 5)
ggplot(data = pd) + geom_point(aes(x = PC1, y = PC2, color = study)) + theme_classic()
dev.off()



