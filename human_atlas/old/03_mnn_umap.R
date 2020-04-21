setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/human_atlas/result/')
data = readRDS('./genebycell.rds')
cellMeta <- readRDS('./cellMeta.rds')
batch = cellMeta$batch
suppressMessages(library(scran))
## highly variable genes
cn <- colnames(data)
fit <- trendVar(data)
decomp <- decomposeVar(data,fit)
diff <- as.numeric(decomp[,'total']) - as.numeric(decomp[,'tech'])
names(diff) = rownames(decomp)
diff = sort(diff[diff>0],decreasing = T)
diff = diff[1: min(1000,length(diff))]
gs <- names(diff)
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
d <- d[,1:50]
library(umap)
set.seed(12345)
u <- umap(d)$layout
saveRDS(u, './mnn_umap.rds')

