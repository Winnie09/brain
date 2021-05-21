library(RColorBrewer)
library(pheatmap)
library(gplots)
library(reshape2)
library(ggplot2)
library(gplots)
dlist.bak <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/data/36nonNormal_combined_log2norm_dlist.rds')
meta <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/data/36nonNormal_combined_meta.rds')
meta.bak <- meta[meta[, 'Pathology'] == 'GBM' & meta[, 'Tumor.Grade'] == 'IV' & meta[, 'Treatment'] == 'Untreated', ]
# for (loc in c( "Temporal", "Parietal", "Frontal")){ ##'Occipital',
######## loc 1
loc1 <- as.character(commandArgs(trailingOnly = T)[[1]])
# loc1 <- 'Frontal'
loc <- loc1
print(loc)
meta <- meta.bak[meta.bak[,'Location'] == loc, ]
allp <- unique(meta[, 'Sample.ID'])
dlist <- dlist.bak[allp]
gene <- lapply(names(dlist), function(p){
  names(which(rowMeans(dlist[[p]]>0) > 0.05))
})
gene2 <- table(unlist(gene))
gene2 <- names(gene2)[gene2==length(dlist)]
gene2 <- gene2[!grepl('^MT-',gene2)]
gene2 <- gene2[!grepl('^RPL|^RPS',gene2)]

for (n in names(dlist)) dlist[[n]] <- dlist[[n]][gene2,]
mat <- do.call(cbind,dlist)
# mat <- mat-rowMeans(mat)  ## center gene expression
mat1 <- (mat-rowMeans(mat))/apply(mat, 1, sd)  ## center gene expression

########### loc 2
loc2 <- as.character(commandArgs(trailingOnly = T)[[2]])
# loc2 <- 'Occipital'
loc <- loc2
print(loc)
meta <- meta.bak[meta.bak[,'Location'] == loc, ]
allp <- unique(meta[, 'Sample.ID'])
dlist <- dlist.bak[allp]
gene <- lapply(names(dlist), function(p){
  names(which(rowMeans(dlist[[p]]>0) > 0.05))
})
gene2 <- table(unlist(gene))
gene2 <- names(gene2)[gene2==length(dlist)]
gene2 <- gene2[!grepl('^MT-',gene2)]
gene2 <- gene2[!grepl('^RPL|^RPS',gene2)]

for (n in names(dlist)) dlist[[n]] <- dlist[[n]][gene2,]
mat <- do.call(cbind,dlist)
# mat <- mat-rowMeans(mat)  ## center gene expression
mat2 <- (mat-rowMeans(mat))/apply(mat, 1, sd)  ## center gene expression

################
int <- intersect(rownames(mat1), rownames(mat2))
str(int)
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/myfunc/01_function.R')
cm <- corfunc(mat1[int, ], mat2[int, ])
saveRDS(cm, paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cellcellcor/res/correlationmatrix_', loc1, '_', loc2, '_centered.rds'))

cellorder1 <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cellcellcor/res/cellorder_', loc1, '.rds')) ##### add
cellorder2 <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cellcellcor/res/cellorder_', loc2, '.rds')) ##### add

## change some outliers values
diag(cm) <- 0
cm[cm > quantile(cm,0.999)] <- quantile(cm,0.999)
p1 <- sub('_.*', '', cellorder1)
p2 <- sub('_.*', '', cellorder2)

png(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cellcellcor/plot/cellcor_', loc1, '_', loc2, '_centered.png'), width = 2000, height = 2000, res = 100)
heatmap(cm[cellorder1, cellorder2], scale = 'none', Rowv = NA, Colv = NA, labRow = NA, labCol = NA,col=colorRampPalette(c('blue3','blue1','skyblue','cyan','green','lightgreen','yellow','orange','pink','red1','red3'))(100),add.expr=eval({abline(v=cumsum(rle(p2)$lengths), lwd=5);abline(h=cumsum(rle(p1)$lengths), lwd=5)}))
dev.off()


