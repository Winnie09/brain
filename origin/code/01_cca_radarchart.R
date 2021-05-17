pdir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/origin/plot/'
rdir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/pseudobulk/data/'
mouseatlaspb <- readRDS(paste0(rdir, 'mouseatlas_pb.rds'))
gbmpb <- readRDS(paste0(rdir, 'gbm_pb.rds'))
gbm <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/integrate/harmony/36nonNormalSeuratGene2000/res/humanAtlas_harmony.rds')
meta <- gbm@meta.data
str(meta)
meta <- meta[meta[,'Pathology'] == 'GBM' & meta[,'Treatment'] == 'Untreated' & meta[,'Tumor.Grade'] == 'IV', ]
meta <- meta[, c('Sample.ID', 'Location')]
meta <- meta[!duplicated(meta), ]

## intersect gene names
int <- intersect(rownames(mouseatlaspb),rownames(gbmpb))
gbmpb <- gbmpb[int,meta[,1]]
mouseatlaspb <- mouseatlaspb[int, ]

mouseatlaspb <- mouseatlaspb[rowSums(mouseatlaspb>1)>1, ]
gbmpb <- gbmpb[rowSums(gbmpb>1)>1, ]

u <- apply(mouseatlaspb, 1, sd)/rowMeans(mouseatlaspb)
v <- apply(gbmpb, 1, sd)/rowMeans(gbmpb)
int <- intersect(names(sort(u[u>median(u)], decreasing = T)), names(sort(v[v>median(v)], decreasing = T)))

gbmpb <- (gbmpb - rowMeans(gbmpb))/apply(gbmpb, 1, sd)
mouseatlaspb <- (mouseatlaspb - rowMeans(mouseatlaspb))/apply(mouseatlaspb, 1, sd)

vmat <- sapply(1:ncol(mouseatlaspb), function(j){
    sapply(1:ncol(gbmpb), function(i){
    cor(gbmpb[int, i], mouseatlaspb[int, j], method = 'spearman')
  })
})
vmat <- as.data.frame(vmat)
dimnames(vmat) <- list(colnames(gbmpb), colnames(mouseatlaspb))

allloc <- meta[match(rownames(vmat),meta[,1]), 2]
pdf(paste0(pdir, 'radarchart_spearman_correlation.pdf'), width = 7, height = 7)
par(mfrow = c(2,2))
for (loc in unique(allloc)){
  tmp <- vmat[allloc == loc, ]
  tmp <- rbind(rep(max(tmp),14), rep(min(tmp),14), tmp)
  radarchart(tmp)  
}
dev.off()

library(pheatmap)
library(gridExtra)
library(ggplot2)
pdf(paste0(pdir, 'hm_spearman_correlation.pdf'), width = 10, height = 10)
plist <- list()
plist[[1]] <- pheatmap(as.matrix(vmat[allloc == sort(unique(allloc))[1], ]))[[4]]
plist[[3]] <- pheatmap(as.matrix(vmat[allloc == sort(unique(allloc))[2], ]))[[4]]
plist[[5]] <- pheatmap(as.matrix(vmat[allloc == sort(unique(allloc))[3], ]))[[4]]
plist[[7]] <- pheatmap(as.matrix(vmat[allloc == sort(unique(allloc))[4], ]))[[4]]
plist[[2]] <- plist[[4]] <- plist[[6]] <- ggplot(data=NULL) + geom_blank() + theme_void()
grid.arrange(grobs = plist,layout_matrix=matrix(c(1,2,3,4,5,6,7),nrow=7))
dev.off()


## CCA coefficient
library(CCA)
int <- intersect(rownames(mouseatlaspb),rownames(gbmpb))
gbmpb <- gbmpb[int,meta[,1]]
mouseatlaspb <- mouseatlaspb[int, ]
res <- cc(mouseatlaspb,gbmpb)
## radar plot
library(fmsb)
# standardize <- function(matrix){ ## standardisze the columns
#   (matrix - colMeans(matrix))/apply(matrix, 2, sd)
# }
# 
# corfunc <- function(m1,m2) {
#   ## calculate the correlation between each column of m1 and m2. 
#       scalematrix(t(m1)) %*% t(scalematrix(t(m2)))/(nrow(m1)-1)      
# }
# 
# scalematrix <- function(data) {
#   cm <- rowMeans(data)
#   csd <- apply(data,1,sd)
#   (data - cm) / csd
# }
# 
u = res$xcoef
v = res$ycoef
# m <- corfunc(standardize(mouseatlaspb[int, ]) %*% u, standardize(gbmpb[int, ]) %*% v)
# m <- as.data.frame(m)
# colnames(m) <- colnames(mouseatlaspb)
# radarchart(m) 

vmat <- as.data.frame(v)
colnames(vmat) <- colnames(mouseatlaspb)
rownames(vmat) <- colnames(gbmpb)

allloc <- meta[match(rownames(vmat),meta[,1]), 2]
par(mfrow = c(2,2))
for (loc in unique(allloc)){
  tmp <- vmat[allloc == loc, ]
  tmp <- rbind(rep(2,14), rep(-2,14), tmp)
  radarchart(tmp)  
}


  
# tests of canonical dimensions
rho <- res$cor
## Define number of observations, number of variables in first set, and number of variables in the second set.
n <- nrow(mouseatlaspb)
p <- ncol(mouseatlaspb)
q <- ncol(gbmpb)

## Calculate p-values using the F-approximations of different test statistics:
library(CCP)
res.ccp <- p.asym(rho, n, p, q, tstat = "Wilks")
p.adjust(res.ccp$p.value, 'fdr')

## correlations wbetween the two sets of variables
mc <- matcor(mouseatlaspb, gbmpb)
mc <- mc[[3]][colnames(gbmpb), colnames(mouseatlaspb)]
mc <- as.data.frame(mc)
allloc <- meta[match(rownames(mc),meta[,1]), 2]
pdf(paste0(pdir, 'radarchart_cca_correlation_samples.pdf'), width = 7, height = 7)
par(mfrow = c(2,2))
for (loc in unique(allloc)){
  tmp <- mc[allloc == loc, ]
  tmp <- rbind(rep(max(mc),14), rep(min(mc),14), tmp)
  mycolor <- colorRampPalette(c('red','cyan','darkblue','purple'))(sum(allloc==loc))
  radarchart(tmp, pcol = mycolor, pfcol=scales::alpha(mycolor,0.1), caxislabels = as.vector(round(quantile(as.vector(unlist(tmp))),2)), cglcol="grey", cglty=1, axislabcol="grey", cglwd=0.8,vlcex=0.8 )
}
dev.off()

## mean across all samples of the same location
mc.mean <- t(sapply(sort(unique(allloc)), function(loc){
  colMeans(mc[allloc == loc, ])
}))
mc.mean <- as.data.frame(mc.mean)
mc.mean <- rbind(rep(max(mc.mean),14), rep(min(mc.mean), 14), mc.mean)
library(RColorBrewer)
mycolor <- colorRampPalette(c('red','purple', 'cyan','darkblue'))(4)
pdf(paste0(pdir, 'radarchart_cca_correlation_location_mean.pdf'), width = 5, height = 5)
radarchart(mc.mean, pcol = mycolor,pfcol=scales::alpha(mycolor,0.1), caxislabels = as.vector(round(quantile(as.vector(unlist(mc.mean))),2)), cglcol="grey", cglty=1, axislabcol="grey", cglwd=0.8,vlcex=0.8 )
dev.off()

