library(Seurat)
library(ggplot2)
library(scattermore)
library(RColorBrewer)
library(ggrepel)
library(here)
# setwd(here())
setwd('/Users/wenpinhou/Dropbox/brain')
cnv = readRDS('atlasGBM/GBMonly/cnvsel/cnv/summary/summary.rds')
alls = sub('_.*', '', names(cnv))


## subset 21 GBM untreated
m <- readRDS('atlasGBM/GBMonly/integrate/seurat/36nonNormalSeuratGene2000/res/meta.rds')
m <- unique(m[,-c(1:4)])
rownames(m) <- m$Sample.ID
m <- m[, c(4:7, 10:18)]
sampsel = rownames(m[m$Pathology == 'GBM' & m$Tumor.Grade == 'IV' & m$Treatment == 'Untreated', ])
str(sampsel)
m = m[sampsel, !colnames(m) %in% c('Pathology', 'Tumor.Grade','Treatment', 'IDH.status', "CDKN2A..2G.loss")]
m$Sex = ifelse(m$Sex == 'F', 1, 0)
m$MGMT.status = ifelse(m$MGMT.status == 'methylated', 1, 0)
m$EGFR.amplification = ifelse(m$EGFR.amplification == 'amplified', 1, 0)
m$PDGFRA.amplification = ifelse(m$PDGFRA.amplification == 'amplified', 1, 0)
m$TP53.mutation = ifelse(m$TP53.mutation == 'mutated', 1, 0)
m$EGFR.mutation = ifelse(m$EGFR.mutation == 'mutated', 1, 0)
m$PTEN.mutation = ifelse(m$PTEN.mutation == 'mutated', 1, 0)

str(alls)
cnv = cnv[alls %in% sampsel]
str(cnv)
alls = sub('_.*', '', names(cnv))

## identify the cnv on the tree
tb = table(cnv, alls)
tb = tb/colSums(tb)
dim(tb)
tb[1:3, 1:3]
tmp = tb > 0.05
summary(rowSums(tmp))

tb2 = tb[which(rowSums(tb > 0.05)>1), ]
str(tb2)


cnv.uni = rownames(tb2)
len = sapply(cnv.uni, function(i){
  length(strsplit(i,';')[[1]])
})
length(cnv.uni)
str(len)
table(len)

trlist = list()
trlist[[1]] = c("chr10-", "chr10-;chr13-", "chr10-;chr13-;chr7+")
trlist[[2]] = c("chr10-", "chr10-;chr22-", "chr10-;chr22-;chr7+")
trlist[[3]] = c("chr10-", "chr10-;chr7+", "chr10-;chr21+;chr7+")
trlist[[4]] = c("chr10-", "chr10-;chr7+", "chr10-;chr22-;chr7+")
trlist[[5]] = c('chr18+')
trlist[[6]] = c('chr7+', "chr10-;chr7+", "chr10-;chr22-;chr7+")
trlist[[7]] = c('chr13-')
trlist[[8]] = c('chr21+')


## celltype
clu <- readRDS('atlasGBM/GBMonly/integrate/seurat/36nonNormalSeuratGene2000//res/cluster.rds')
libsize = table(sub('_.*', '', names(clu)))
ct <- read.table('doc/GBM_cluster_annotations.csv', as.is = TRUE, sep = ',', header = T)
ct <- ct[, c(1,3)]
ct <- ct[match(clu, as.character(ct[,1])),2]
ct.gbm <- ct
names(ct.gbm) <- names(clu)
str(ct.gbm)
str(cnv)
ct.gbm = ct.gbm[names(cnv)]
df = data.frame(cell = names(ct.gbm), sample = sub('_.*','', names(ct.gbm)), ct = ct.gbm, cnv = cnv,  stringsAsFactors = FALSE)
rownames(df) = df[,1]


## for each sample, the cell proportion of each of 13 cnv 

s = unique(df$sample)[1]

df.s = df[df[,2] == s, ]
df.s = df[df[,4] %in% trlist[[1]], ]
table(df.s[,2], df.s[])


branchprop = matrix(0, nrow = length(sampsel), ncol = length(trlist))
dimnames(branchprop) = list(sampsel, paste0('branch', 1:length(trlist)))

for (i in 1:length(trlist)){
  tr.tmp = trlist[[i]]
  v = table(df[df[,4] %in% tr.tmp, 2])
  branchprop[names(v), i] = v/libsize[names(v)]
}

pdir = 'atlasGBM/GBMonly/cnvtree/plot/'
pdf(paste0(pdir, 'hm_branchprop_on_branch_by_sample.pdf'), width = 5, height = 6)
pheatmap(branchprop)
dev.off()


## for each cnv pattern, test sample covariate and ct association
m.tmp = m
cellprop = branchprop
pval = matrix(NA, nrow = ncol(m.tmp), ncol = ncol(cellprop))
dimnames(pval) = list(colnames(m.tmp), colnames(cellprop))

for (tmpct in colnames(pval)){
  for (tmpmeta in rownames(pval)){
    print(paste0(tmpct, ';', tmpmeta))
    fstat = summary(lm(cellprop[, tmpct] ~ model.matrix(~m.tmp[, tmpmeta])[,-1]))$fstatistic
    pval[tmpmeta, tmpct] = pf(fstat[1],fstat[2],fstat[3],lower.tail = F) ## tail
  }
}

fdrmat = matrix(p.adjust(as.vector(pval)), nrow = nrow(pval), ncol = ncol(pval))
dimnames(fdrmat) = dimnames(pval)

## plot FDR heatmap
mat_breaks <- c(0, 0.001, 0.01, 0.05, 0.1, 0.3,0.5,1)
library(viridis)
pdf(paste0(pdir, 'fdr_branch_by_meta.pdf'), width = 5, height = 4.2)
pheatmap(fdrmat, color = inferno(length(mat_breaks) - 1), breaks = mat_breaks, show_rownames = T, show_colnames = T)
dev.off()

## plot barplots
for (id in c(2,3,4,8)){
pd = rbind(data.frame(prop = branchprop[m$TP53.mutation == '1', paste0('branch',id)], type = 'mutated'),
           data.frame(prop = branchprop[m$TP53.mutation == '0', paste0('branch',id)], type = 'unmutated'))
pdf(paste0(pdir, paste0('branch',id), '_TP53_boxplot.pdf'), width = 3, height = 2.1)
print(ggplot(data = pd, aes(x = type, y = prop, color = type)) +
  geom_boxplot() +
  geom_point(size = 0.3)+
  theme_classic() + xlab('TP53') + ylab(paste0('branch',id, ' proportion')) +
  theme(legend.position = 'none'))
dev.off()  
}

