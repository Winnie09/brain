result_path = as.character(commandArgs(trailingOnly = T)[[1]])
key = as.character(commandArgs(trailingOnly = T)[[2]])
print(key)
print('result_path')
if(!exists('result_path')) print('please provide a working directory path')
setwd(result_path)
# result_path ='/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlas/res/'
# key = 'sub'
# key = 'full'
# result_path = '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/human_atlas/result/'
# key = 'referenceVec'
# key = 'referenceVec_allenDEG'

library(Seurat)
library(ggplot2)
library(cowplot)
library(pheatmap)
library(gplots)
library(RColorBrewer)

meta = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/GBM/old/Nov9/res/seurat/meta.rds')
res = readRDS(paste0('./originCt/result/',key,'/originCt_prop.rds'))
m = t(res[,grepl('GBM',colnames(res))])
colnames(m) = paste0('clu',res[,'cluster'],'(',res[,'atlas_celltype'],')')
meta[,5] = as.factor(meta[,5])
meta[,6] = as.factor(meta[,6])
meta[,7] = as.factor(meta[,7])
meta1 = meta[match(sub('\\(.*','',rownames(m)),sub('\\(.*','',meta$sample)), 5:7]
row.names(meta1) <- row.names(m)
meta2 = read.csv('/home-4/whou10@jhu.edu/scratch/Wenpin/GBM/doc/Meta_Wenpin_10292019.csv')  
meta2$Sample.ID = sapply(meta2$Sample.ID, function(i) sub(' ','',i))
meta2 = meta2[match(sub('\\(.*','',rownames(m)), meta2$Sample.ID),]
meta2 = meta2[,!colnames(meta2)%in%c('Progression.date')]

metadata = cbind(meta1, meta2[,2:16])
rownames(meta1) <- rownames(meta2) <- row.names(metadata) <- row.names(m)
rownames(res) <- colnames(m)
set.seed(12345)
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)
morecols2 <- colorRampPalette(brewer.pal(6,"Set3"))
breaksList = seq(0, 1, by = 0.05)

my_color = sapply(colnames(meta2)[2:ncol(meta2)], function(i) {
  if (is.factor(meta2[,i])){
    n = length(levels(meta2[,i]))
    col = morecols2(n)
    names(col) = levels(meta2[,i])
  } else {
    n = length(unique(meta2[,i])) 
    col = morecols2(n)
    names(col) = unique(meta2[,i])
  }
  col  
})
dir.create(paste0('./originCt/plot/',key), showWarnings = F, recursive = T)
pdf(paste0('./originCt/plot/',key,'/hm_batch.pdf'),width=10,height=10)
dir.create(paste0('./originCt/plot/',key),showWarnings = F, recursive = T)
pheatmap(m, show_rownames=T, cluster_cols=T, cluster_rows=T, scale="row",
         cex=1, border_color=FALSE,annotation_row=meta1,
         color = rev(morecols(20)),breaks=breaksList,
         annotation_col=res[,c('gbm_enrichment','atlas')])
dev.off()

pdf(paste0('./originCt/plot/',key,'/hm_clinical.pdf'),width=12,height=8)
pheatmap(m, show_rownames=T, cluster_cols=T, cluster_rows=T, scale="row",
         annotation_colors=my_color[2:14],
         cex=1, border_color=FALSE,annotation_row=meta2[,2:9],
         color = rev(morecols(20)),breaks=breaksList)
dev.off()


pdf(paste0('./originCt/plot/',key,'/hm_clinical2.pdf'),width=12,height=8)
pheatmap(m, show_rownames=T, cluster_cols=T, cluster_rows=T, scale="row",
         annotation_colors=my_color[2:14],
         cex=1, border_color=FALSE,annotation_row=meta2[,10:16],
         color = rev(morecols(20)),breaks=breaksList)
dev.off()

