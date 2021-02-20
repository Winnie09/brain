library(here)
setwd(here())
pdir <- 'atlasGBM/mouseGBM/integrate/harmony/36nonNormal_seuratGene/plot/'

library(Seurat)
res <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/mouseGBM/integrate/harmony/36nonNormal_seuratGene/res/atlasGBM_harmony_with_ct.rds')
data = res@assays$RNA@data
harmony <- res@reductions$harmony@cell.embeddings
clu <- as.numeric(res@active.ident) ## clu 0 will be 1 on the plot
names(clu) <- rownames(harmony)
ct <- res@meta.data$celltype
names(ct) <- rownames(res@meta.data)
ct.clu <- unlist(sapply(min(clu):max(clu), function(i){
  tab <- table(ct[names(clu)[which(clu == i)]])
  n <- names(tab)[which.max(tab)[1]]
  if (length(n) > 0){
    n
  } else{
    'NA'
  }
}))
names(ct.clu) <- min(clu):max(clu)
library(TSCAN)
library(igraph)
library(RColorBrewer)
set.seed(12345)
mc <- exprmclust(t(harmony),cluster=clu,reduce=F)
# ord <- TSCANorder(mc,orderonly = T)
png(paste0(pdir, 'MST_vertexLabeledAsCelltype_only.png'), width = 2000, height = 2000, res = 200)
set.seed(12345)
V(mc$MSTtree)$color <- brewer.pal(9, 'Blues')[2]
V(mc$MSTtree)$size <- 1
plot(mc$MSTtree, vertex.label = ct.clu, vertex.label.cex = 0.8,
     vertex.label.color="black", edge.color = 'black', edge.width=2,
     label.dist = 1)
dev.off()

png(paste0(pdir, 'MST_vertexLabeledAsCluster.png'), width = 2000, height = 2000, res = 200)
set.seed(12345)
V(mc$MSTtree)$color <- brewer.pal(9, 'Blues')[2]
V(mc$MSTtree)$size <- 6
plot(mc$MSTtree, vertex.label.color="black", edge.color = 'black', edge.width=2,label.dist = 1)
dev.off()


# ord <- TSCANorder(mc,orderonly = T)


plotMST <- function(mst, gene, geneByCellExpr, clu, vertex.label = clu){
  ## mst: an igraph object. For example, mc$MSTtree where mc is the result of TSCAN::exprmclust(mat).
  ## clu: a vector or length = num.cells. cell clusters.
  expr <- geneByCellExpr[gene, ]
  clumean <- tapply(expr, clu, mean)
  V(mst)$color <- brewer.pal(9, 'Blues')[5]
  V(mst)$size <- clumean * 5
  set.seed(12345)
  plot(mst, vertex.label = vertex.label, vertex.label.color="black", edge.color = 'black', edge.width=2, main = gene)
}


tb = read.table('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/doc/Marker_Genes.csv',header=T,sep=',',as.is=T)
library(gridExtra)
for (i in seq(1,ncol(tb))){
  print(i)
  plist = list()
  gvec = intersect(tb[,i], rownames(data))
  for (g in gvec){
    print(g)
    png(paste0(pdir,paste0('MST_',colnames(tb)[i],'_markerGene_',g,'.png')), width = 1800, height = 1800, res = 200)
    print(plotMST(mc$MSTtree, g, data, clu, vertex.label = ct.clu))
    dev.off()
  }
}




