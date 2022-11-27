rm(list=ls())
#setwd('/home/whou10/scratch16/whou10/brain/')
setwd('/Users/wenpinhou/Dropbox/brain/')
## suva's data, grouped by the chetan data cell clusters 
ddir <- 'predDNAm/integrate/mapquery/'
ddir.ref <- 'atlasGBM/GBMonly/integrate/seurat/36nonNormalSeuratGene2000/res/'
u <- readRDS(paste0(ddir, 'umap.rds'))
colnames(u) <- c('x', 'y')
u.ref <- readRDS(paste0(ddir.ref, 'umap_embeddings.rds'))
ct <- readRDS(paste0(ddir.ref, 'celltype.rds'))


di <- spatstat.geom::crossdist(u[,1], u[,2], u.ref[,1], u.ref[,2])
dimnames(di) <- list(rownames(u), rownames(u.ref))
str(di)

ctassign <- sapply(rownames(di), function(c){
  tab <- table(ct[names(sort(di[c, ])[1:10])])
  names(which.max(tab))
}) ## alternatively, can fit a multivariate normal
table(ctassign)
str(ctassign)
saveRDS(ctassign, 'predDNAm/knn/ctassign/suva_data_ctassign.rds')

## =============================================
## plot chetan's and suva's data using ctassign
## =============================================
pd <- data.frame(ct = ctassign, u1 = u[names(ctassign), 1], u2 = u[names(ctassign), 2])
str(pd)
library(ggplot2)
pdf('predDNAm/knn/ctassign/ctassign_text_repel.pdf', width =30, height = 30)
ggplot(data = pd, aes(x = u1, y = u2, label = ct, color = ct))+
  geom_point() + 
  ggrepel::geom_text_repel( ) + 
  theme_classic()
dev.off()


## =======================
## predict and evaluation 
## =======================
##### in each cluster, randomly pick one. use others' mean to predict this one. 
## compared predicted and measure profiles: PCC, SCC, MSE 
ddir <- 'predDNAm/data/sep/GSE151506_hg38/proc/scRRBS/'
d1 <- readRDS(paste0(ddir, 'promoter_up1k_down1k_filtered_AC_MES.rds'))
d2 <- readRDS(paste0(ddir, 'promoter_up1k_down1k_filtered_NPC_OPC.rds'))
int = intersect(rownames(d1), rownames(d2))
d1 <- d1[int, ]
d2 <- d2[int, ]
d <- cbind(d1, d2)
d.sp <- sub('_.*', '', colnames(d))
pred.list <- true.list <- list()
for (i in unique(ctassign)) {
  print(i)
    gp <- intersect(names(ctassign)[ctassign == i], colnames(d))
    if (length(gp) < 5) next
    sp <- sub('_.*','', gp)
    pred <- as.data.frame(do.call(cbind, lapply(unique(sp), function(s){
      testid <- gp[sp == s] ## leave one sample out
      trainid <- setdiff(gp, testid)
      pred.tmp <- rowMeans(d[, trainid, drop=FALSE], na.rm = TRUE)  
    })))
    true <- as.data.frame(do.call(cbind, lapply(unique(sp), function(s){
      testid <- gp[sp == s] ## leave one sample out
      true.tmp <- rowMeans(d[, testid, drop=FALSE], na.rm = TRUE)
    })))
    colnames(true) <- colnames(pred) <- paste0(unique(sp), ';', i)
    true.list[[i]] <- true
    pred.list[[i]] <- pred
}
names(true.list) <- names(pred.list) <- NULL
pred <- do.call(cbind, pred.list)
true <- do.call(cbind, true.list)


##### randomness
pred.rd.list  <- list()
for (i in unique(ctassign)) {
    print(i)
    gp <- intersect(names(ctassign)[ctassign == i], colnames(d))
    if (length(gp) < 5) next
    sp <- sub('_.*','', gp)
    pred.tmp <- as.data.frame(do.call(cbind, lapply(unique(sp), function(s){
      testid <- gp[sp == s] ## leave one sample out
      set.seed(123)
      trainid <- sample(setdiff(colnames(d)[d.sp != s], gp), length(gp) - length(testid))
      rowMeans(d[, trainid, drop=FALSE], na.rm = TRUE)  
    })))
    colnames(pred.tmp) <- paste0(unique(sp), ';', i)
    pred.rd.list[[i]] <- pred.tmp
}
names(pred.rd.list) <- NULL
pred.rd <- do.call(cbind, pred.rd.list)
identical(colnames(pred), colnames(pred.rd))
## evaluation functions
scalematrix <- function(data) {
  cm <- rowMeans(data)
  csd <- apply(data,1,sd)
  (data - cm) / csd
}

corfunc <- function(m1,m2) {
  ## calculate the correlation between each column of m1 and m2. 
  scalematrix(t(m1)) %*% t(scalematrix(t(m2)))/(nrow(m1)-1)      
}

mycor<- function(x, method, na.action){
  ## calculate the correlation between columns of x, ignoring NAs.
  r<- apply(x, 2, function(j){
    apply(x, 2, function(i){
      as.numeric(cor.test(i,j, method = method, na.action = na.action)$estimate)
    })
  })
  out<-c()
  out$r<- r
  return(out) 
}

## evaluation
alls <- sub(';.*', '', colnames(pred))
table(alls)

LOO <- c()
for (s in sort(unique(alls))){
  print(s)
  pred.s <- pred[, alls == s] ###############
  true.s <- true[, alls == s]
  
  id <- which(complete.cases(cbind(true.s, pred.s))) 
  if (length(id) == 0) next
  cross.sp <- corfunc(t(true.s[id,]), t(pred.s[id,]))
  cross.sp <- cross.sp[1,] ## median 0.25
  # cross.sp <- sapply(rownames(pred.s), function(g){
  #   id <- which(!is.na(pred.s[g, ]) & !is.na(true.s[g, ]))
  #   cor(as.vector(as.matrix(pred.s[g, id])), as.vector(as.matrix(true.s[g, id])))
  # })  ## median 0.03
  
  # install.packages('corbetw2mat')
  
  cross.pro<- mycor(cbind(true.s, pred.s) , method="pearson", na.action=na.omit)
  cross.pro <- cross.pro[[1]][1,]
  pred.rd.s <- pred.rd[, alls == s] ###################
  id <- which(complete.cases(cbind(true.s, pred.rd.s))) ## 0
  cross.sp.rd<- corfunc(t(true.s[id,]), t(pred.rd.s[id,]))
  cross.sp.rd <- cross.sp.rd[1,]
  
  # cross.sp.rd <- summary(sapply(rownames(pred.rd.s), function(g){
  #   id <- which(!is.na(pred.rd.s[g, ]) & !is.na(true.s[g, ]))
  #   cor(as.vector(as.matrix(pred.rd.s[g, id])), as.vector(as.matrix(true.s[g, id])))
  # }))  ## median 0.09
  cross.pro.rd <- mycor(cbind(true.s, pred.rd.s) , method="pearson", na.action=na.omit)
  cross.pro.rd <- cross.pro.rd[[1]][1,]
  # str(cross.sp)
  # str(cross.sp.rd)
  # 
  # str(cross.pro)
  # str(cross.pro.rd)
  # 
  # summary(cross.sp)
  # summary(cross.sp.rd)
  # 
  # summary(cross.sp)
  # summary(cross.sp.rd[names(cross.sp)])
  # 
  # summary(cross.pro)
  # summary(cross.pro.rd)
  

  cross.sp.var <- sapply(1:nrow(true.s), function(i){
    tmp <- true.s[i, ]
    var(tmp[!is.na(tmp)])
  })
  names(cross.sp.var) <- rownames(true.s)
  LOO[[s]] <- list(cross.sp = cross.sp, cross.pro = cross.pro, 
                   cross.sp.rd = cross.sp.rd, cross.pro.rd = cross.pro.rd,
                   cross.sp.var = cross.sp.var)
}
saveRDS(LOO, 'predDNAm/knn/cv/res/LOO_pcc.rds')

str(LOO)

### =============
### plot results
### =============
pdir <- 'predDNAm/knn/cv/plot/'
source('/Users/wenpinhou/Dropbox/resource/ggplot_theme.R')
#source('/home/whou10/scratch16/whou10/resource/ggplot_theme.R')
theme_set(.new_theme)
library(RColorBrewer)



## cross sample




pd <- lapply(names(LOO), function(s){
  cross.sp.var <- LOO[[s]][[5]]
  names(cross.sp.var) <- rownames(true.s)
  a <- quantile(cross.sp.var, na.rm = T)
  pro.sel <- names(cross.sp.var)[cross.sp.var >= a[4]]
  tmp <- rbind(data.frame(sample = paste0(s, '(',table(alls)[s],')'), type = 'data', cross.sp = LOO[[s]][[1]]), 
               data.frame(sample = paste0(s, '(',table(alls)[s],')'), type = 'random', cross.sp = LOO[[s]][[3]]))
  int <- intersect(rownames(tmp), pro.sel)
  tmp <- tmp[int, ]
  
})
str(pd)
pd <- do.call(rbind, pd)

pdf(paste0(pdir, 'LOO_cross.sp_4th_quantile.pdf'), width = 5, height = 1.5)
ggplot(data = pd, aes(x = sample, y = cross.sp, color = type, fill = type)) +
  geom_boxplot(alpha = 0.2) + 
  #geom_violin(alpha  = 0.2, width = 0.2, scale = 'width') +
  #geom_point(size = 0.2, stroke = 0, alpha = 0.2) + 
  scale_color_brewer(palette = 'Set1') + 
  xlab('Leave-out Sample') + 
  ylab('Cross-sample PCC')
dev.off()


## cross promoter
pd <- lapply(names(LOO), function(s){
  tmp <- rbind(data.frame(sample = paste0(s, '(',table(alls)[s],')'), type = 'data', cross.pro = LOO[[s]][[2]]), 
               data.frame(sample = paste0(s, '(',table(alls)[s],')'), type = 'random', cross.pro = LOO[[s]][[4]]))
  
})
str(pd)
pd <- do.call(rbind, pd)

pdf(paste0(pdir, 'LOO_cross_promoter_at_least_5cell_in_ct.pdf'), width = 5, height = 1.5)
ggplot(data = pd, aes(x = sample, y = cross.pro, color = type, fill = type)) +
  geom_boxplot(alpha = 0.2) + 
  #geom_violin(alpha  = 0.2, width = 0.2, scale = 'width') +
  #geom_point(size = 0.2, stroke = 0, alpha = 0.2) + 
  scale_color_brewer(palette = 'Set1') + 
  xlab('Leave-out Sample') + 
  ylab('Cross-sample PCC')
dev.off()

