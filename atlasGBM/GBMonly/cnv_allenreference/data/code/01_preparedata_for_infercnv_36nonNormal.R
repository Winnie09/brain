library(Seurat)
library(here)
setwd(here())

## prepare reference: select allen brain 10x cells 
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/myfunc/01_function.R')
atlas <-  countfunc(readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/allen_human_10x/proc/count.rds'), platform='10x', log2normalize = TRUE, species='human')

colnames(atlas) <- gsub('allen:', 'allen_human_10x:', colnames(atlas))

meta <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/humanatlas/data/meta_allcell.rds')
meta[,1] <- sub('allen_human_10x;', '', meta[,1])
rownames(meta) <- meta[,1]
meta <- meta[colnames(atlas), ]

selcell <- lapply(unique(meta[,7]), function(ct){
  tmp <- rownames(meta[which(meta[,7] == ct), ])
  if (length(tmp) > 500){
    tmp[1:500]
  } else {
    tmp
  }
})  
selcell <- unlist(selcell)  
atlas <- atlas[, selcell] ##
meta <- data.frame(cell = selcell, cluster = meta[selcell,7], stringsAsFactors = F) ##

## prepare query: select GBM cells
brain <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/integrate/harmony/36nonNormalSeuratGene2000/res/humanAtlas_harmony.rds')
clu <- paste0('cluster ', as.character(Idents(brain)))
names(clu) <- colnames(brain@assays$RNA@counts)
# clu[clu == 'cluster 5'] <- 'oligodendrocyte'
# clu[clu == 'cluster 9'] <- 'microglia'
# clu[clu == 'cluster 10'] <- 'endothelial'
# clu[clu == 'cluster 13'] <- 'oligodendrocyte'
# clu[clu == 'cluster 15'] <- 'oligodendrocyte'
patient = brain@meta.data[,'Sample.ID']
rm(brain)

expr <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/data/36nonNormal_combined_count_mat.rds')
int <- intersect(rownames(expr), rownames(atlas))
str(int)
expr <- expr[int, ]
atlas <- atlas[int, ]


## gene annotation file
library(data.table)
gtf <- fread('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/gtf/grch38.gtf',data.table=F)  ## grch38 = hg38（13161, more genes intersect with hg19）
gtf <- gtf[gtf[,3]=='gene',]
gn <- sub('\";.*','',sub('.*; gene_name \"','',gtf[,9]))
gr <- data.frame(gene = row.names(expr),gtf[match(row.names(expr),gn),c(1,4,5)],stringsAsFactors = F)
colnames(gr) = c('gene','chr','start','end')
gr <- gr[!is.na(gr[,4]),]
write.table(gr,file='atlasGBM/GBMonly/cnv_allenreference/data/data/36nonNormal_allen10x_combined_gr.txt',quote=F,sep='\t',col.names=F,row.names=F)


p = unique(patient)[1]
for (p in unique(patient)){
  print(p)
  id.p <- which(patient == p)
  tab <- table(clu[id.p])
  id <- which(clu %in% names(tab)[tab > 2])
  id <- intersect(id, id.p)
  ## cell annotation file
  anno <- data.frame(cell = names(clu)[id], cluster =  clu[id], stringsAsFactors = FALSE)
  anno <- rbind(anno, meta)
  write.table(anno, paste0('atlasGBM/GBMonly/cnv_allenreference/data/data/36nonNormal_allen10x_combined_', p, '_cellanno.txt'), quote=F, sep="\t", row.names = FALSE, col.names = FALSE)  ##
  ## write gene by cell expression matrix
  k <- expr[, id]
  k <- k[gr[,1], ]
  df = as.data.frame(cbind(k, atlas[rownames(k), ]), stringsAsFactors=F)
  fwrite(df, file = paste0('atlasGBM/GBMonly/cnv_allenreference/data/data/36nonNormal_allen10x_combined_', p, '_log2NormMatrix.txt'), quote = F, sep = '\t', row.names = T)
}
