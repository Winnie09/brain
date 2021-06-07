library(NMF)
library(parallel)

dlist <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/data/36nonNormal_combined_log2norm_dlist.rds')
meta <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/data/36nonNormal_combined_meta.rds')
meta <- meta[meta[, 'Pathology'] == 'GBM' & meta[, 'Tumor.Grade'] == 'IV' & meta[, 'Treatment'] == 'Untreated', ]
path <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/nmf/res/overallgenepathway.rds')

# ==============
# gene filtering 
# ==============
allp <- unique(meta[, 'Sample.ID'])
gene <- lapply(allp, function(p){
  mattmp <- dlist[[p]]          
  rownames(mattmp[rowMeans(mattmp>0) > 0.01, ])
})
gene2 <- table(unlist(gene))
gene2 <- names(gene2)[gene2==length(allp)]
gene2 <- gene2[!grepl('^MT-',gene2)]
gene2 <- gene2[!grepl('^RPL|^RPS',gene2)]

ave <- rowMeans(sapply(allp, function(i){
  rowMeans(dlist[[i]][gene2,])
}))
genegroup <- cut(ave, quantile(ave,seq(0,1,length.out=26)),include.lowest = T)
names(genegroup) <- names(ave)

centerfunc <- function(a) a-rowMeans(a)

factorscore <- sapply(1:ncol(path),function(i) {
  gs <- path[,i]
  gcont <- as.vector(sapply(gs, function(g){
    group <- genegroup[g]
    set.seed(12345)
    sample(setdiff(names(genegroup[genegroup == group]),gs), 100, replace = F)
  }))
  unlist(lapply(allp,function(n) colMeans(dlist[[n]][gs,]) - colMeans(dlist[[n]][gcont,])))
})
saveRDS(factorscore, '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/nmf/res/pathwayscore.rds')


