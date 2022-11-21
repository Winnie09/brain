
## suva's data, grouped by the chetan data cell clusters 
ddir <- '/home/whou10/scratch16/whou10/brain/predDNAm/integrate/mapquery/'
ddir.ref <- '/home/whou10/scratch16/whou10/brain/atlasGBM/GBMonly/integrate/seurat/36nonNormalSeuratGene2000/res/'
u <- readRDS(paste0(ddir, 'umap.rds'))
colnames(u) <- c('x', 'y')
u.ref <- readRDS(paste0(ddir.ref, 'umap_embeddings.rds'))
ct <- readRDS(paste0(ddir.ref, 'celltype.rds'))

tmp1 <- data.frame(x = tapply(u.ref[,1], ct, mean),
                   y = tapply(u.ref[,2], ct, mean))
tmp2 <- data.frame(x = tapply(u.ref[,1], ct, max),
                   y = tapply(u.ref[,2], ct, max))
tmp3 <- data.frame(x = tapply(u.ref[,1], ct, min),
                   y = tapply(u.ref[,2], ct, min))
tmp1 <- cbind(celltype = rownames(tmp1), tmp1)
tmp2 <- cbind(celltype = rownames(tmp2), tmp2)
tmp3 <- cbind(celltype = rownames(tmp2), tmp3)
center <- rbind(tmp1, tmp2, tmp3)

d <- rbind(u, center[, 2:3])
di <- as.matrix(dist(d))
ctassign <- sapply(rownames(u), function(i){
  tmp <- t(t(center[,2:3]) - u[i,])
  center[which.min(rowSums(tmp^2)), 1]
})
str(ctassign)
table(ctassign)

## in each cluster, randomly pick one. use others' mean to predict this one. 
## compared predicted and measure profiles: PCC, SCC, MSE 
ddir <- '/home/whou10/scratch16/whou10/brain/predDNAm/data/sep/GSE151506_hg38/proc/scRRBS/'
d1 <- readRDS(paste0(ddir, 'promoter_up1k_down1k_filtered_AC_MES.rds'))
d2 <- readRDS(paste0(ddir, 'promoter_up1k_down1k_filtered_NPC_OPC.rds'))
int = intersect(rownames(d1), rownames(d2))
d1 <- d1[int, ]
d2 <- d2[int, ]
d <- cbind(d1, d2)
str(d)

perf <- as.data.frame(do.call(rbind, lapply(unique(ctassign), function(i){
    print(i)
    gp <- names(ctassign)[ctassign == i]
    gp <- intersect(gp, colnames(d))
    set.seed(123)
    testid.all <- sample(gp, min(10, length(gp)))
    tmp <- t(sapply(testid.all, function(testid){
      trainid <- setdiff(gp, testid)
      pred <- rowMeans(d[, trainid, drop=FALSE], na.rm = TRUE)
      id <- which(!is.na(d[,testid, drop=FALSE]) & !is.na(pred))
      c(cor(d[id,testid], pred[id]), cor(d[id,testid], pred[id], method = 'spearman'), sum((d[id,testid] - pred[id])^2))
    }))
    colnames(tmp) <- c('PCC', 'SCC', 'MSE')
    tmp <- cbind(celltype = i, tmp)
})))
perf[,2] <- as.numeric(perf[,2])
perf[,3] <- as.numeric(perf[,3])
perf[,4] <- as.numeric(perf[,4])
str(perf)


## in random case. still select those testing samples. 
## randomly select training samples 
perf.rd <- as.data.frame(do.call(rbind, lapply(unique(ctassign), function(i){
    print(i)
    gp <- names(ctassign)[ctassign == i]
    gp <- intersect(gp, colnames(d))
    set.seed(1234)
    testid.all <- sample(gp, min(10, length(gp)))
    tmp <- t(sapply(testid.all, function(testid){
      set.seed(1234)
      trainid <- sample(colnames(d),  length(gp) - length(testid))
      pred <- rowMeans(d[, trainid, drop=FALSE], na.rm = TRUE)
      id <- which(!is.na(d[,testid, drop=FALSE]) & !is.na(pred))
      c(cor(d[id,testid], pred[id]), cor(d[id,testid], pred[id], method = 'spearman'), sum((d[id,testid] - pred[id])^2))
    }))
    colnames(tmp) <- c('PCC', 'SCC', 'MSE')
    tmp <- cbind(celltype = i, tmp)
})))
perf.rd[,2] <- as.numeric(perf.rd[,2])
perf.rd[,3] <- as.numeric(perf.rd[,3])
perf.rd[,4] <- as.numeric(perf.rd[,4])
str(perf.rd)
summary(perf.rd[,2])
summary(perf[,2])
summary(perf.rd[,3])
summary(perf[,3])
summary(perf.rd[,4])
summary(perf[,4])


