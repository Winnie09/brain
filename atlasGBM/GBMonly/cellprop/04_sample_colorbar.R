rm(list=ls())
library(Seurat)
library(ggplot2)
library(scattermore)
library(RColorBrewer)
library(pheatmap)
library(here)
# setwd(here())
setwd('/Users/wenpinhou/Dropbox/brain')
tb = readRDS('atlasGBM/GBMonly/cellprop/cellprop_by_samples.rds')

m <- readRDS('atlasGBM/GBMonly/integrate/seurat/36nonNormalSeuratGene2000/res/meta.rds')
m <- unique(m[,-c(1:4)])
rownames(m) <- m$Sample.ID
m <- m[, c(4:7, 10:18)]

cs = readRDS('atlasGBM/GBMonly/cellprop/color_scheme.rds')
colv = cs$colv
annotation_row = cs$annotation_row
annotation_colors = cs$annotation_colors

head(rownames(m))
m$Pathology

v = rownames(m)
for (i in 1:nrow(m)){
  print(i)
  if (grepl('Oli',m$Pathology[i])) {
    v[i] = sub('GBM', 'OLI', v[i])
  } else if (grepl('Ast',m$Pathology[i])) {
    v[i] = sub('GBM', 'AST', v[i])
  } else {
    next
  }
}

df = data.frame(old.name = rownames(m), new.name = v, stringsAsFactors = F)
rownames(df) = df[,1]
saveRDS(df, 'atlasGBM/GBMonly/cellprop/GBM_samples_rename.rds')
write.csv(df, 'atlasGBM/GBMonly/cellprop/GBM_samples_rename.csv')


rownames(tb) <- df[rownames(tb), 2]
rownames(m) <- df[rownames(m), 2]
str(tb)
str(m)
m = m[order(m[,3],m[,5],m[,2]), ]
tb = tb[rownames(m),]
m = m[, c(3,4, 5, 2, 1, 6:13)]
pdf('atlasGBM/GBMonly/cellprop/sample_colorbar.pdf', width = 10, height = 11)
pheatmap(mat = tb[,1:2],
         color = colv,
         annotation_row = m, 
         annotation_colors = annotation_colors[c(4,2,1,3, 5:11)],
         cluster_rows = F,
         cluster_cols = T, 
         border_color = NA)
dev.off()

