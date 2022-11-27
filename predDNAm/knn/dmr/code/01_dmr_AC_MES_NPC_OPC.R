ddir <- '/home/whou10/scratch16/whou10/brain/predDNAm/data/sep/GSE151506_hg38/proc/scRRBS/'

## retain promoters with > 5 CpG in >= 10 cells per group (NPC/OPC-like, AC/MES-like) 
d <- readRDS(paste0(ddir, 'paired.rds'))
d.bak = d
sample <- sub('_.*', '', colnames(d))
meta.sample <- readRDS('/home/whou10/scratch16/whou10/brain/predDNAm/data/sep/GSE151506_hg38/meta/meta.rds')
d <- d[,sample %in% meta.sample[meta.sample[,3]=='GBM', 1]]
ct <- sub('\\..*', '', sub(':.*', '', sub('.*_', '', colnames(d))))

library(GenomicRanges)
library(data.table)
g <- fread('/home/whou10/scratch16/whou10/resource/gtf/grch38.gtf',data.table=F)
g <- g[g[,3]=='gene',]
gn <- sub('".*','',sub('.*gene_name "','',g[,9]))
ddir <- rdir <- '/home/whou10/scratch16/whou10/brain/predDNAm/data/sep/GSE151506_hg38/proc/scRRBS/'
cr <- GRanges(seqnames=g[,1],IRanges(start=g[,4],end=g[,5]),strand=g[,7])
names(cr) <- gn
pro <- promoters(cr,upstream=1000,downstream=1000)

for (i in 1:2){
  if (i == 1){
    ct.sel <- c('AC', 'MES')
    ct.save <- 'AC_MES'
  } else if (i == 2){
    ct.sel <- c('NPC', 'OPC')
    ct.save <- 'NPC_OPC'
  }
  print(ct.sel)
  m <- d[, ct %in% ct.sel]
  l <- as.numeric(sub('.*_','',rownames(m)))
  mr <- GRanges(seqnames=sub('_.*','',rownames(m)),IRanges(start=l,end=l))
  o <- as.matrix(findOverlaps(mr,pro))
  sp <- split(o[,1], gn[o[,2]])
  s <- sapply(sp,function(i) {
    if (length(i) > 5){
      tmp <- colMeans(m[i,,drop=F],na.rm=T)
      if (sum(!is.na(tmp)) > 10){
        return(tmp)
      } else NA
    } else NA
  })
  s <- do.call(rbind, s[!is.na(s)])
  print(str(s))
  saveRDS(s,file=paste0(rdir, '/promoter_up1k_down1k_filtered_', ct.save, '.rds'))
}

## mean DNAme across all regions
s1 <- readRDS(file=paste0(rdir, '/promoter_up1k_down1k_filtered_AC_MES.rds'))
s2 <- readRDS(file=paste0(rdir, '/promoter_up1k_down1k_filtered_NPC_OPC.rds'))
int = intersect(rownames(s1), rownames(s2))
s1 <- s1[int, ]
s2 <- s2[int, ]
global.mean1 <- colMeans(s1, na.rm = T)
global.mean2 <- colMeans(s2, na.rm = T)
  
  
## lm of a promoter DNAme ~ group + mean DNAme
## P-values were calculated from the t-statistic (computed by dividing the cellular state coefficient by the standard error).
pval <- sapply(int, function(i){
  df <- rbind(data.frame(me = s1[i, ], ct = 'AC_MES', global.mean = global.mean1),
  data.frame(me = s2[i, ], ct = 'NPC_OPC', global.mean = global.mean2))
  mod <- lm(me~ct+global.mean, data = df)
  summary(mod)$coefficients[2,4]
})
sum(pval < 0.05, na.rm = T)
fdr <- p.adjust(pval)
sum(fdr < 0.05, na.rm = T)

## volcano plot
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


