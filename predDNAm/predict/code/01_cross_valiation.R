## calculate standard score matrix 
scalematrix <- function(data) {
  cm <- rowMeans(data)
  csd <- sqrt((rowMeans(data*data) - cm^2) / (ncol(data) - 1) * ncol(data))
  (data - cm) / csd
}

## calculate pearson correlation between columns
corfunc <- function(m1,m2,type='concordant') { 
  if (type=='concordant') {
    rowSums(scalematrix(m1) * scalematrix(m2),na.rm=T)/(ncol(m1)-1)
  } else {
    scalematrix(t(m1)) %*% t(scalematrix(t(m2)))/(nrow(m1)-1)            
  }
}

source('/home/whou10/scratch16/whou10/software/metpred/predict.R')
source('/home/whou10/scratch16/whou10/software/metpred/trainmodel.R')
pdir <- '/home/whou10/scratch16/whou10/brain/predDNAm/predict/plot/'
rdir <- '/home/whou10/scratch16/whou10/brain/predDNAm/predict/res/'
d <- readRDS('/home/whou10/scratch16/whou10/brain/predDNAm/combat_pca/res/combat_parametric_adjustment.rds')
m <- readRDS('/home/whou10/scratch16/whou10/brain/predDNAm/data/combine/res/combinemecgi.rds')

d <- d[,colnames(m)]
m <- m[rowMeans(is.na(m))==0,]
stu <- sub(':.*','',colnames(d))
names(stu) <- colnames(d)

## cross-valiation
k <- sapply(unique(stu),function(i) {
  testn <- names(which(stu==i))
  trainn <- names(which(stu!=i))
  mod <- trainmodel(d[,trainn],m[,trainn])  
  pred <- predict(d[,testn],mod)
  tm <- m[,testn]
  correlation <- corfunc(pred,tm) # correlation: standarize 
  csd <- apply(tm,1,sd)
  cbind(correlation,csd)
},simplify = F) 
saveRDS(k, paste0(rdir, 'cross_validation_correlation_columnSd.rds'))

