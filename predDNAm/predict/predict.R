setwd('/home/whou10/scratch16/whou10/brain/predDNAm/data/filter')
source('/hpc/group/jilab/zj/metgbm/software/trainmodel.R')
source('/hpc/group/jilab/zj/metgbm/software/predict.R')
g <- readRDS('ge.rds')
m <- readRDS('me.rds')
mod1 <- trainmodel(g$geo,m$geo)
pred1 <- predict(g$tcga,mod1)
cv1_row <- sapply(1:nrow(pred1),function(i) cor(m$tcga[i,],pred1[i,]))
pdf('/home/whou10/scratch16/whou10/brain/predDNAm/predict/geotrain_acrosssample.pdf')
hist(cv1_row)
dev.off()

cv1_col <- sapply(1:ncol(pred1),function(i) cor(m$tcga[,i],pred1[,i]))
pdf('/home/whou10/scratch16/whou10/brain/predDNAm/predict/geotrain_acrosssite.pdf')
hist(cv1_col)
dev.off()

# summary(cv1)
# sd1 <- apply(m$tcga,1,sd)
# summary(cv1[sd1 > 0.1])

mod2 <- trainmodel(g$tcga,m$tcga)
pred2 <- predict(g$geo,mod2)

cv2_row <- sapply(1:nrow(pred2),function(i) cor(m$geo[i,],pred2[i,]))
pdf('/home/whou10/scratch16/whou10/brain/predDNAm/predict/tcgatrain_acrosssample.pdf')
hist(cv2_row)
dev.off()

cv2_col <- sapply(1:ncol(pred2),function(i) cor(m$geo[,i],pred2[,i]))
pdf('/home/whou10/scratch16/whou10/brain/predDNAm/predict/tcgatrain_acrosssite.pdf')
hist(cv2_col)
dev.off()



