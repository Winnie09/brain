#setwd('/home/whou10/scratch16/whou10/brain/')
setwd('/Users/wenpinhou/Dropbox/brain/')
ddir <- rdir <- 'predDNAm/data/sep/GSE151506_hg38/proc/scRRBS/'

## read in data
s1 <- readRDS(file=paste0(rdir, '/promoter_up1k_down1k_filtered_AC_MES.rds'))
s2 <- readRDS(file=paste0(rdir, '/promoter_up1k_down1k_filtered_NPC_OPC.rds'))
int = intersect(rownames(s1), rownames(s2))
s1 <- s1[int, ]
s2 <- s2[int, ]
d <- cbind(s1,s2)

## ctassign
ctassign <- readRDS('predDNAm/knn/ctassign/suva_data_ctassign.rds')
ctassign <- ctassign[colnames(d)]
tab <- table(ctassign)
ctassign <- ctassign[ctassign %in% names(tab)[tab > 1]]
table(ctassign)
d <- d[, names(ctassign)]
str(d)

## mean DNAme across all regions
dlist <- global.mean <- list()
for (i in unique(ctassign)){
  print(i)
  dlist[[i]] <- d[, ctassign == i, drop = F]
  global.mean[[i]] <- colMeans(dlist[[i]], na.rm = T)
}
str(dlist)
str(global.mean)

  
## make sure each group in differential has enough (e.g. >=5) non-NA cells
## lm of a promoter DNAme ~ group + mean DNAme
## P-values were calculated from the t-statistic (computed by dividing the cellular state coefficient by the standard error).
pvallist <- list()
for (ct1id in 1 : (length(unique(ctassign)) - 1)){
  ct1 <- unique(ctassign)[ct1id]
  for (ct2id in (ct1id + 1) : (length(unique(ctassign)))){
    ct2 <- unique(ctassign)[ct2id]
    print(paste0(ct1, ';', ct2))
    pvallist[[paste0(ct1, ';', ct2)]] <- sapply(int, function(i){
      me1 <- dlist[[ct1]][i, ]
      me2 <- dlist[[ct2]][i, ]
      if (sum(!is.na(me1)) < 5 | sum(!is.na(me2)) < 5) return(NA)
      me1 <- me1[!is.na(me1)]
      me2 <- me2[!is.na(me2)]
      df <- rbind(data.frame(me = me1, ct = 1, global.mean = global.mean[[ct1]][names(me1)]),
                  data.frame(me = me2, ct = 0, global.mean = global.mean[[ct2]][names(me2)]))
      df <- df[complete.cases(df), ,drop=F]
      mod <- lm(me~ct+global.mean, data = df)
      # if (nrow(summary(mod)$coefficients) == 1) {
      #   return(NA)
      # } else {
        summary(mod)$coefficients[2,4]
      # }
    })  
  }
}

saveRDS(pvallist, 'predDNAm/knn/dmr/res/suva_ctassign_pairwise_dmr_pval.rds')

ctcomp <- data.frame(ct1 = sub(';.*','', names(pvallist)), ct2 = sub('.*;','', names(pvallist)))

pval <- do.call(cbind, pvallist)
dim(pval)
fdr <- apply(pval, 2, p.adjust)
# fdr <- lapply(colnames(pval), function(i){
#   tmp <- pval[,i]
#   tmp <- tmp[!is.na(tmp)]
#   summary(p.adjust(tmp))
#   summary(p.adjust(pval[,i]))
# })
str(fdr)

numdiff <- sapply(colnames(fdr), function(i){
  sum(fdr[,i] < 0.05, na.rm = T)
})
summary(numdiff)
str(numdiff)
table(numdiff)

dmr <- sapply(colnames(fdr), function(i){
  rownames(fdr)[which(fdr[,i] < 0.05)]
})
str(dmr)
dmr <- dmr[sapply(dmr, length)>0]
saveRDS(dmr, 'predDNAm/knn/dmr/res/suva_ctassign_pairwise_dmr_promoter.rds')

dmr.m <- sapply(names(dmr), function(i){
  if (length(dmr[[i]]) < max(sapply(dmr, length))) {
    c(dmr[[i]], rep(' ', max(sapply(dmr, length)) - length(dmr[[i]])))
  } else {
    dmr[[i]]
  }
})
write.csv(dmr.m, 'predDNAm/knn/dmr/res/suva_ctassign_pairwise_dmr_promoter.csv')


## plot number of DMR for each ct pair



## gene ontology of the dmr in each pair of cts comparison


## volcano plot for each pair of cts' dmr 
diff <- rowMeans(s1, na.rm = T) - rowMeans(s2, na.rm = T)
pd <- data.frame(mediff = diff, p = -log10(pval), stringsAsFactors = F)
g <- c('CDX2', 'FGF3', 'HOXD8', 'FGF5', 'GATA6', 'POU4F2', 'FOXL1', 'ESPN', 'GABRA4', 'FOXD2')
pd2 <- cbind(pd[g, ],label = g)
pd.diff <- pd[pd[,2] > -log10(0.05) & abs(pd[,1])>0.05, ]

library(ggplot2)
source('/home/whou10/scratch16/whou10/resource/ggplot_theme.R')
theme_set(.new_theme)
pdir <- '/home/whou10/scratch16/whou10/brain/predDNAm/knn/dmr/plot/'
pdf(paste0(pdir, 'volcano.pdf'), width = 2.5, height = 2.5)
ggplot() + 
  geom_point(data = pd, aes(x = mediff, y = p), alpha = 0.05, size = 0.1, )+
  geom_point(data = pd.diff, aes(x = mediff, y = p), size = 0.5, alpha = 0.3)+
  geom_point(data = pd2, aes(x = mediff, y = p), color = 'red', size = 0.8)+
  ggrepel::geom_text_repel(data = pd2, aes(x = mediff, y = p, label = label), size = 2, color = 'red')+
  geom_vline(xintercept=c(-0.05, 0.05), col="black", linetype = 'dashed') +
  geom_hline(yintercept=-log10(0.05), col="black", linetype = 'dashed') +
  xlab('DNA methylation difference') + 
  ylab(expression(-log[10](P)))+
  ggtitle(paste0(length(int), ' promoters (+/-1k bp), L:AC/MES hypo'))
dev.off()



