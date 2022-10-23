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
stu <- sub(':.*', '', colnames(combine))

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

## ------------------
## pca for each study
## ------------------
### --------------------------------------------------------------------------
for (i in c('JHU', 'GSE119776', 'GSE121723')){
  print(i)
  rdir <- '/home/whou10/scratch16/whou10/brain/predDNAm/pca/'
  dir.create(paste0(rdir, i), recursive = T)
  PCA(genebycellmat = combine[, stu==i], save.pca = T, plot.statistics=F, 
      plot.dir = paste0(rdir, i, '/'),
      result.dir = paste0(rdir, i, '/'),
      PC_for_columns = TRUE,
      findVariableGenes = TRUE, findVariableGenesMethod = 'lm', maxVariableGenes = 2e3, numPC = 10, smoothFittingMethod = 'loess')
  pr  <- readRDS(paste0(rdir, i, '/pr.rds'))
  str(pr)
  prfull <- readRDS(paste0(rdir, '/prfull.rds'))
  str(prfull)
  pcvar = (prfull[[1]])^2
  pcvar <- pcvar/sum(pcvar)
  sum(pcvar)
  
  alls <- sub('.*:', '', rownames(pr))
  head(alls)
  table(alls)
  pd = data.frame(PC1 = pr[,1], PC2 = pr[,2], stringsAsFactors = F, sample = alls)
  rownames(pd) <- rownames(pr)
  table(pd$sample)
  library(ggplot2)
  library(ggrepel)
  library(RColorBrewer)
  source('/home/whou10/scratch16/whou10/resource/ggplot_theme.R')
  theme_set(.new_theme) 
  
  pdf(paste0(rdir, i, '/pca_sample.pdf'), height = 3, width = 3)
  print(ggplot(
    data = pd,
    aes(x = PC1, y = PC2, color = sample, label = sample)
  ) +
    geom_point(size = 0.6) +
    scale_color_manual(values = colorRampPalette(rev(
      brewer.pal(8, 'Set1')
    ))(length(unique(
      alls
    )))) +
    geom_text_repel(
      size = 2,
      max.overlaps = Inf,
      segment.size = 0.05
    ) + 
    theme(legend.position = 'none') + 
    xlab(paste0('Principal component 1 (', round(pcvar[1],4)*100, '%)')) +
    ylab(paste0('Principal component 2 (', round(pcvar[2],4)*100, '%)')))
  dev.off()
  
}

### --------------------------------------------------------------------------
i = 'GSE151506'
dir.create(paste0(rdir, i), recursive = T)
PCA(genebycellmat = combine[, stu==i], save.pca = T, plot.statistics=F, 
    plot.dir = paste0(rdir, i, '/'),
    result.dir = paste0(rdir, i, '/'),
    PC_for_columns = TRUE,
    findVariableGenes = TRUE, 
    findVariableGenesMethod = 'lm', 
    maxVariableGenes = 2e3, numPC = 10, 
    smoothFittingMethod = 'loess')
pr  <- readRDS(paste0(rdir, i, '/pr.rds'))
prfull <- readRDS(paste0(rdir, '/prfull.rds'))
pcvar = (prfull[[1]])^2
pcvar <- pcvar/sum(pcvar)
alls <- gsub('_.*', '', sub('.*:', '', rownames(pr)))
allct <- gsub('.*_', '', sub('.*:', '', rownames(pr)))
pd = data.frame(PC1 = pr[,1], PC2 = pr[,2],
                sample = alls, 
                sample.label = alls,
                celltype = allct,
                celltype.label = allct,
                stringsAsFactors = F)
rownames(pd) <- rownames(pr)
pd[duplicated(pd[,4]),4] <- NA
pd[duplicated(pd[,4]),6] <- NA
for (colvar in c('sample', 'celltype')){
  pdf(paste0(rdir, i, '/pca_', colvar, '.pdf'), height = 3, width = 3)
  assign('var', colvar)
  assign('var.label', paste0(colvar, '.label'))
  print(ggplot(
    data = pd,
    aes_string(x = 'PC1', y = 'PC2', color = var, label = var.label)
  ) +
    geom_point(size = 0.6) +
    scale_color_manual(values = colorRampPalette(rev(
      brewer.pal(8, 'Set1')
    ))(length(unique(
      pd[,var]
    )))) +
    geom_text_repel(
      size = 2,
      max.overlaps = Inf,
      segment.size = 0.05
    ) + 
    theme(legend.position = 'none') + 
    xlab(paste0('Principal component 1 (', round(pcvar[1],4)*100, '%)')) +
    ylab(paste0('Principal component 2 (', round(pcvar[2],4)*100, '%)')))
  dev.off()
}

### --------------------------------------------------------------------------
table(stu)
i = 'GSE100351'
dir.create(paste0(rdir, i), recursive = T)
PCA(genebycellmat = combine[, stu==i], save.pca = T, plot.statistics=F, 
    plot.dir = paste0(rdir, i, '/'),
    result.dir = paste0(rdir, i, '/'),
    PC_for_columns = TRUE,
    findVariableGenes = TRUE, 
    findVariableGenesMethod = 'lm', 
    maxVariableGenes = 2e3, numPC = 10, 
    smoothFittingMethod = 'loess')
pr  <- readRDS(paste0(rdir, i, '/pr.rds'))
str(pr)
prfull <- readRDS(paste0(rdir, '/prfull.rds'))
pcvar = (prfull[[1]])^2
pcvar <- pcvar/sum(pcvar)

alls <- gsub('_', '', substr(sub('.*:', '', rownames(pr)), 9, 12))
names(alls) <- rownames(pr)
head(alls)
length(unique(alls))
meta <- readRDS('/home/whou10/scratch16/whou10/brain/predDNAm/data/sep/GSE100351_hg38/meta/meta.rds')
str(meta)
alls2 <- sub('.*:', '', rownames(pr))
str(alls2)
length(intersect(alls2, meta$Description))
pd = data.frame(PC1 = pr[,1], PC2 = pr[,2],
                sample = alls, 
                sample.label = alls,
                stringsAsFactors = F)
rownames(pd) <- rownames(pr)
pd[duplicated(pd[,4]),4] <- NA
for (colvar in c('sample')){
  pdf(paste0(rdir, i, '/pca_', colvar, '.pdf'), height = 3, width = 3)
  assign('var', colvar)
  assign('var.label', paste0(colvar, '.label'))
  print(ggplot(
    data = pd,
    aes_string(x = 'PC1', y = 'PC2', color = var, label = var.label)
  ) +
    geom_point(size = 0.6) +
    scale_color_manual(values = colorRampPalette(rev(
      brewer.pal(8, 'Set1')
    ))(length(unique(
      pd[,var]
    )))) +
    geom_text_repel(
      size = 2,
      max.overlaps = Inf,
      segment.size = 0.05
    ) + 
    theme(legend.position = 'none') + 
    xlab(paste0('Principal component 1 (', round(pcvar[1],4)*100, '%)')) +
    ylab(paste0('Principal component 2 (', round(pcvar[2],4)*100, '%)')))
  dev.off()
}

### --------------------------------------------------------------------------
table(stu)
i = 'TCGA'
dir.create(paste0(rdir, i), recursive = T)
PCA(genebycellmat = combine[, stu==i], save.pca = T, plot.statistics=F, 
    plot.dir = paste0(rdir, i, '/'),
    result.dir = paste0(rdir, i, '/'),
    PC_for_columns = TRUE,
    findVariableGenes = TRUE, 
    findVariableGenesMethod = 'lm', 
    maxVariableGenes = 2e3, numPC = NULL, 
    smoothFittingMethod = 'loess')
pr  <- readRDS(paste0(rdir, i, '/pr.rds'))
str(pr)
prfull <- readRDS(paste0(rdir, '/prfull.rds'))
pcvar = (prfull[[1]])^2
pcvar <- pcvar/sum(pcvar)
sum(pcvar)
alls <- gsub('.*-', '', sub('.*:', '', rownames(pr)))
names(alls) <- rownames(pr)
head(alls)
length(unique(alls))

pd = data.frame(PC1 = pr[,1], PC2 = pr[,2],
                sample = alls, 
                sample.label = alls,
                stringsAsFactors = F)
rownames(pd) <- rownames(pr)
pd[duplicated(pd[,4]),4] <- NA
for (colvar in c('sample')){
  pdf(paste0(rdir, i, '/pca_', colvar, '.pdf'), height = 3, width = 3)
  assign('var', colvar)
  assign('var.label', paste0(colvar, '.label'))
  print(ggplot(
    data = pd,
    aes_string(x = 'PC1', y = 'PC2', color = var, label = var.label)
  ) +
    geom_point(size = 0.6) +
    scale_color_manual(values = colorRampPalette(rev(
      brewer.pal(8, 'Set1')
    ))(length(unique(
      pd[,var]
    )))) +
    geom_text_repel(
      size = 2,
      max.overlaps = Inf,
      segment.size = 0.05
    ) + 
    theme(legend.position = 'none') + 
    xlab(paste0('Principal component 1 (', round(pcvar[1],4)*100, '%)')) +
    ylab(paste0('Principal component 2 (', round(pcvar[2],4)*100, '%)')))
  dev.off()
}


