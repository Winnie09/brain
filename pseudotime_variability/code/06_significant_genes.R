setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/pseudotime_variability/result')
allfn = list.files('./6_10_3_pm')
lenv = NULL
for (fn in allfn){
  lenv[fn] = length(list.files(paste0('./6_10_3_pm/', fn)))
}
lenv = lenv[lenv==21]
lenv = lenv[!names(lenv)%in%c(61,62,70,73,80,87)]
pval = 0
for (fn in names(lenv)){
  print(fn)
  p = readRDS(paste0('./6_10_3_pm/', fn,'/pvalue.rds'))
  pval = p + pval
}
pval = pval/100/length(lenv)
for (i in 1:ncol(pval)){
  pval[, i] = p.adjust(pval[, i], method='fdr')
}
sig <- (pval<0.05) + 0
genes <- rownames(sig)[which(rowSums(sig) > 0)]
df = data.frame(gene = rownames(sig),num.significant.base = rowSums(sig,na.rm = T), stringsAsFactors = F)
df = df[order(df[,2],decreasing = T),]
saveRDS(df,'./6_10_3/significant_genes.rds')


