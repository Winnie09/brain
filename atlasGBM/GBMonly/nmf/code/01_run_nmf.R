library(NMF)
library(parallel)
num.nmf <- 4
rdir <- paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/nmf/res/', num.nmf, '/')
dir.create(rdir, showWarnings = F, recursive = T)

dlist <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/data/36nonNormal_combined_log2norm_dlist.rds')
meta <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/data/36nonNormal_combined_meta.rds')
meta <- meta[meta[, 'Pathology'] == 'GBM' & meta[, 'Tumor.Grade'] == 'IV' & meta[, 'Treatment'] == 'Untreated', ]
# i = unique(meta[, 'Location'])[1]
# allp <- unique(meta[meta[,'Location']==i, 'Sample.ID'])

# ==========================
# only use maglinant cells !!
# =========================



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

# ==================================
# center expr for each gene, and nmf
# ==================================
res <- mclapply(allp,function(p) {
  mat <- dlist[[p]][gene2,]
  mat <- mat - rowMeans(mat) ## centered expression
  mat[mat < 0]  <-  0  ## convert negatives to zero
  nmf(mat, num.nmf, seed = 12345) ## nmf 
},mc.cores = 4)
names(res) <- allp
saveRDS(res, paste0(rdir, 'res.rds'))

### ====================
### get gene signatures
### ====================
num.nmf <- 4
rdir <- paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/nmf/res/', num.nmf, '/')

allres = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/nmf/res/res.rds')
dlist <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/data/36nonNormal_combined_log2norm_dlist.rds')
dlist <- dlist[names(allres)]
library(NMF)
## get 
allscore <- do.call(cbind,lapply(names(allres), function(p){
  res = allres[[p]]
  res <- basis(res)
  colnames(res) <- paste0(p, ';factor', 1:ncol(res))
  res
}))


ave <- rowMeans(sapply(dlist, function(i){
  rowMeans(i)
}))
genegroup <- cut(ave, quantile(ave,seq(0,1,length.out=26)),include.lowest = T)
names(genegroup) <- names(ave)

factorscore <- apply(allscore,2,function(i) {
  gs <- names(sort(i, decreasing = T)[1:30])
  
  gcont <- as.vector(sapply(gs, function(g){
    group <- genegroup[g]
    set.seed(12345)
    sample(names(genegroup[genegroup == group]), 100, replace = F)
  }))
  unlist(lapply(names(allres),function(n) colMeans(dlist[[n]][gs,]) - colMeans(dlist[[n]][gcont,])))
  
})
saveRDS(factorscore, '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/nmf/res/factorscore.rds')


## cluster factor pathways
library(pheatmap)
hclu <- cutree(hclust(as.dist(1-cor(factorscore))),k = 4)

overallgenepathway <- sapply(1:max(hclu),function(i) {
  tmpclu <- names(hclu)[hclu==i]
  
  names(sort(rowMeans(allscore[,tmpclu,drop=F]), decreasing = T)[1:30])
}) 
saveRDS(overallgenepathway, '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/nmf/res/overallgenepathway.rds')

factorscore = factorscore[, names(sort(hclu))]


library(RColorBrewer)
png('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/nmf/plot/factorscore_select.png', width = 1500, height = 1500, res = 100)
mycolor <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100)
mycolor <- c(rep(mycolor[1],30), rep(mycolor[1:30],each = 10), mycolor[31:69], rep(mycolor[70:99], each = 10), rep(mycolor[length(mycolor)], 30))
pheatmap(t(factorscore), scale = 'row', show_rownames = T, show_colnames = F, cluster_cols = F, cluster_rows = F, color = mycolor)
dev.off()


