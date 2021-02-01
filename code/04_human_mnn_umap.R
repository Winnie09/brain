data = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/result/merge_human/genebycell.rds')
batch <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/result/merge_human/batch.rds')
batch = rep(names(batch), batch)
suppressMessages(library(scran))
## highly variable genes
cn <- colnames(data)
fit <- trendVar(data)
decomp <- decomposeVar(data,fit)
gs <- row.names(decomp)[decomp[,'total'] > decomp[,'tech']]
data <- data[gs,]
scmd <- sapply(1:length(unique(batch)),function(sp) {
  paste0("data[gs,batch==unique(batch)[",sp,"]]")
})
# run mnn
d <- min(50,nrow(data))
ncores <- 4
cmd <- paste0("res <- fastMNN(",paste0(scmd,collapse = ","),",BPPARAM=MulticoreParam(ncores),d=",d,')')
eval(parse(text=cmd))
d <- res$corrected
row.names(d) <- cn
# d <- d[,1:30]
library(umap)
set.seed(12345)
u <- umap(d)$layout
saveRDS(u, '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/result/merge_human/mnn_umap.rds')
