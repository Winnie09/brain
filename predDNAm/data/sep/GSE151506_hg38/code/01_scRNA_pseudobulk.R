rm(list=ls())
ddir <- '/home/whou10/scratch16/whou10/brain/predDNAm/data/sep/GSE151506_hg38/raw/scRNA/'
rdir <- '/home/whou10/scratch16/whou10/brain/predDNAm/data/sep/GSE151506_hg38/proc/scRNA/'
af = list.files(ddir)
alld <- lapply(af, function(f){
  print(f)
  d <- read.table(paste0(ddir, f), sep = '\t')
  colnames(d) <- d[1,]
  d <- d[-1,]
  rownames(d) <- d[,1]
  d <- d[,-1]
  dn <- dimnames(d)
  d2 <- matrix(as.numeric(as.matrix(d)), nrow = nrow(d))
  dimnames(d2) <- dn
  d2
})

alld2 <- do.call(cbind, alld)
alld2.bak = alld2
str(alld2)
cn <- colnames(alld2)
for (i in 1:length(cn)){
  tmp <- cn[i]
  cn[i] <- paste0(strsplit(tmp, '\\.')[[1]], collapse = '_')
  cn[i] <- sub('Pl1', 'P1', cn[i])
  cn[i] <- sub('Pl2', 'P2', cn[i])
  cn[i] <- sub('P_1', 'P1', cn[i])
  cn[i] <- sub('P_2', 'P2', cn[i])
}
colnames(alld2) <- cn
saveRDS(alld2, paste0(rdir, 'tpm_allsamples.rds'))

## meta data
mdir <- '/scratch16/hji7/whou10/brain/predDNAm/data/sep/GSE151506_hg38/meta/'
tb1 <- read.csv(paste0(mdir, 'GSE151506_suppltb5_GBM_malignant_classification.csv'))
tb2 <- read.csv(paste0(mdir, 'GSE151506_suppltb5_IDHMutant_malignant_classification.csv'))
selcell <- c(tb1$CellAssignment, tb2$CellAssignment)
names(selcell) <- c(tb1[,1], tb2[,1])
selcell <- selcell[nchar(selcell) > 0]
names(selcell) <- sub('P_1', 'P1', names(selcell))
names(selcell) <- sub('P_2', 'P2', names(selcell))
names(selcell) <- sub('\\.', '-', names(selcell))

## pseudobulk
d <- alld2[, intersect(names(selcell), colnames(alld2))]
selcell <- selcell[colnames(d)]
allp <- sub('_.*', '', colnames(d))
allp <- substr(allp, 1, 6)

pblist <- lapply(sort(unique(allp)), function(p){
  tmp <- d[, allp == p]
  tmp.cell <- selcell[intersect(colnames(tmp), names(selcell))]
  pb <- sapply(sort(unique(tmp.cell)), function(i){
    tmp <- rowMeans(tmp[,names(tmp.cell)[tmp.cell == i], drop = FALSE])
  })
  colnames(pb) <- paste0(p,';',colnames(pb))
  pb
})
pb <- do.call(cbind, pblist)
pb <- log2(pb +1)
saveRDS(pb, paste0(rdir, 'pseudobulk_ge.rds'))
