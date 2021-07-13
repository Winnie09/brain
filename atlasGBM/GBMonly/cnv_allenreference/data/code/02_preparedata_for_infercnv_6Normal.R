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
dlist = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/data/6Normal_combined_log2norm_dlist.rds')
int <- intersect(rownames(dlist[[1]]), rownames(atlas))
str(int)
atlas <- atlas[int, ]


## gene annotation file
library(data.table)
gtf <- fread('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/gtf/grch38.gtf',data.table=F)  ## grch38 = hg38（13161, more genes intersect with hg19）
gtf <- gtf[gtf[,3]=='gene',]
gn <- sub('\";.*','',sub('.*; gene_name \"','',gtf[,9]))
gr <- data.frame(gene = int,gtf[match(int,gn),c(1,4,5)],stringsAsFactors = F)
colnames(gr) = c('gene','chr','start','end')
gr <- gr[!is.na(gr[,4]),]
write.table(gr,file='atlasGBM/GBMonly/cnv_allenreference/data/data/6Normal_allen10x_combined_gr.txt',quote=F,sep='\t',col.names=F,row.names=F)


# p = unique(patient)[1]
# for (p in unique(patient)){
for (p in names(dlist)){
  print(p)
  # id.p <- which(patient == p) ## cell of this sample
  # tab <- table(clu[id.p]) ## cluster of this sample
  # id <- which(clu %in% names(tab)[tab > 2]) ## more than one cell?
  # id <- intersect(id, id.p)
  
  ## cell annotation file
  # anno <- data.frame(cell = names(clu)[id], cluster =  clu[id], stringsAsFactors = FALSE)
  anno <- data.frame(cell = colnames(dlist[[p]]), cluster = 'non-tumor', stringsAsFactors = FALSE)
  anno <- rbind(anno, meta)
  write.table(anno, paste0('atlasGBM/GBMonly/cnv_allenreference/data/data/6Normal_allen10x_combined_', p, '_cellanno.txt'), quote=F, sep="\t", row.names = FALSE, col.names = FALSE)  ##
  ## write gene by cell expression matrix
  k <- dlist[[p]]
  k <- k[gr[,1], ]
  df = as.data.frame(cbind(k, atlas[rownames(k), ]), stringsAsFactors=F)
  fwrite(df, file = paste0('atlasGBM/GBMonly/cnv_allenreference/data/data/6Normal_allen10x_combined_', p, '_log2NormMatrix.txt'), quote = F, sep = '\t', row.names = T)
}

