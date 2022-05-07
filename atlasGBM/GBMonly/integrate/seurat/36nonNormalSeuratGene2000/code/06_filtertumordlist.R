library(Matrix)

dlist <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/integrate/seurat/36nonNormalSeuratGene2000/res/norm.rds')

clu <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/integrate/seurat/36nonNormalSeuratGene2000/res/cluster.rds')
normalclu <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_myelymreference/res/normalclu.rds')
dlist <- dlist[,!clu %in% normalclu]

samp <- sub('_.*','',colnames(dlist))
meta <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/data/36nonNormal_combined_meta.rds')
meta <- meta[meta[, 'Pathology'] == 'GBM' & meta[, 'Tumor.Grade'] == 'IV' & meta[, 'Treatment'] == 'Untreated', ]
allp <- unique(meta[, 'Sample.ID'])
dlist <- sapply(allp,function(s) {
  dlist[,samp==s]
})

gene <- lapply(allp, function(p){
  mattmp <- dlist[[p]]          
  rownames(mattmp[rowMeans(mattmp>0) > 0.01,])
})
gene2 <- table(unlist(gene))
gene2 <- names(gene2)[gene2==length(allp)]
gene2 <- gene2[!grepl('^MT-',gene2)]
gene2 <- gene2[!grepl('^RPL|^RPS',gene2)]
dlist <- sapply(dlist,function(i) i[gene2,])
saveRDS(dlist,file='/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/integrate/seurat/36nonNormalSeuratGene2000/res/filtertumordlist.rds')
