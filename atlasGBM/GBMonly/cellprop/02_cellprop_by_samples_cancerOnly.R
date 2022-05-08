library(Seurat)
library(ggplot2)
library(scattermore)
library(RColorBrewer)
library(pheatmap)
library(here)
# setwd(here())
setwd('/Users/wenpinhou/Dropbox/brain')
clu <- readRDS('atlasGBM/GBMonly/integrate/seurat/36nonNormalSeuratGene2000//res/cluster.rds')
ct <- read.table('doc/GBM_cluster_annotations.csv', as.is = TRUE, sep = ',', header = T)
ct <- ct[, c(1,3)]
ct <- ct[match(clu, as.character(ct[,1])),2]
ct.gbm <- ct
names(ct.gbm) <- names(clu)
ct.gbm = ct.gbm[!ct.gbm %in% c('Microglia', 'CBV', 'Endothelial', 'OLGs')]
df = data.frame(cell = names(ct.gbm), sample = sub('_.*','', names(ct.gbm)), ct = ct.gbm, stringsAsFactors = FALSE)

tb = table(df[,2], df[,3])
str(tb)
tb <- tb/rowSums(tb)
rowSums(tb)
saveRDS(tb, 'atlasGBM/GBMonly/cellprop/cellprop_by_samples_cancerOnly.rds')


m <- readRDS('atlasGBM/GBMonly/integrate/seurat/36nonNormalSeuratGene2000/res/meta.rds')
m <- unique(m[,-c(1:4)])
rownames(m) <- m$Sample.ID
m <- m[, c(4:7, 10:18)]

cs = readRDS('atlasGBM/GBMonly/cellprop/color_scheme.rds')
colv = cs$colv
m = cs$annotation_row
annotation_colors = cs$annotation_colors

mybreak = c(0, seq(0.001, max(tb), length.out = 98))

## ===== change sample names GBM to OLI and AST
df = readRDS('atlasGBM/GBMonly/cellprop/GBM_samples_rename.rds')
rownames(tb) <- df[rownames(tb), 2]
str(tb)


pdf('atlasGBM/GBMonly/cellprop/cellprop_by_samples_cancerOnly.pdf', width = 10, height = 11)
pheatmap(mat = tb,
         color = colv,
         annotation_row = m, 
         annotation_colors = annotation_colors,
         cluster_cols = T, 
         border_color = NA,
         breaks = mybreak)
dev.off()


### two sample t-test

for (i in 1:ncol(m)){
  v = m[,i]
  v[is.na(v)] = 'unknown'
  m[,i] = v
}
pval <- matrix(NA, nrow = ncol(m), ncol = ncol(tb))
dimnames(pval) = list(colnames(m), colnames(tb))
for (tmpct in colnames(pval)){
  for (tmpmeta in rownames(pval)){
    print(paste0(tmpct,'_',tmpmeta))
    fstat <- summary(lm(tb[, tmpct]~model.matrix(~m[, tmpmeta])[,-1]))$fstatistic
    pval[tmpmeta, tmpct] = pf(fstat[1],fstat[2],fstat[3],lower.tail = F) ## tail distribution
  }
}
str(pval)
fdrmat = matrix(p.adjust(as.vector(pval)), nrow = nrow(pval), ncol = ncol(pval))
dimnames(fdrmat) = dimnames(pval)

colv = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
colv = rev(c(rep(colv[1], 20), 
             rep(colv[30], 20), 
             rep(colv[50], 20), 
             rep(colv[71:75], each = 4), 
             rep(colv[96:100], each = 4)))
mybreak = c(seq(0, 0.001, length.out = 20), 
            seq(0.001,0.01, length.out = 21)[-1],
            seq(0.01,0.05, length.out = 21)[-1],
            seq(0.05,0.1, length.out = 21)[-1],
            seq(0.1, 1, length.out = 21)[-1])

pdf('atlasGBM/GBMonly/cellprop/celltype_meta_association_cancerOnly.pdf', width = 5, height = 4.2)
pheatmap(fdrmat, color = colv, breaks = mybreak, show_rownames = T, show_colnames = T)
dev.off()

