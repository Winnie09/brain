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

### ----------------------------------
### plot sample cluster-proportion on MSTtree
clu.tmp <- clu[grepl('GBM', names(clu))]
sp <- sub('_.*', '', names(clu.tmp))
prop <- sapply(seq(min(clu), max(clu)), function(i){
  tmp <- clu.tmp[clu.tmp==i]
  tab <- table(sub('_.*', '', names(tmp)))
  v <- rep(0, length(unique(sp)))
  names(v) <- unique(sp)
  v[names(tab)] <- tab
  v  
})
colnames(prop) <- seq(min(clu), max(clu))
prop <- prop/rowSums(prop)



meta <- res@meta.data
meta <- meta[grepl('GBM', meta$study), 6:31]
meta <- meta[,-c(3:10)]
colnames(meta)
meta <- meta[!duplicated(meta$study), ]
rownames(meta) <- meta$study

for (fea in setdiff(colnames(meta), c('study', 'age', 'MDSC','Mutation.Rate '))){
  print(fea)
  len <- length(unique(meta[,fea]))
  if (len > 4){
    png(paste0(pdir, 'MST_GBM_feature_',fea, '.png'), width = 900*sqrt(len), height = 900*(sqrt(len) + 1), res = 100)
    par(mfrow = c(sqrt(len), (sqrt(len)+1)))
  } else {
    png(paste0(pdir, 'MST_GBM_feature_',fea, '.png'), width = 900*len, height = 900, res = 100)
    par(mfrow = c(1,len))
  }
  unifea <- unique(meta[, fea])
  unifea <- unifea[!is.na(unifea)]
  unifea <- setdiff(unifea,'')
  for (i in unifea){
    id <- which(meta[,fea] == i)
    a = prop[meta[id, 'study'], ,drop=FALSE]
    med <- apply(a, 2, mean)
    mst = mc$MSTtree
    V(mst)$color <- brewer.pal(9, 'Blues')[5]
    V(mst)$size <- med * 100
    set.seed(12345)
    print(plot(mst, vertex.label = ct.clu, vertex.label.color="black", edge.color = 'black', edge.width=2, main = i, label.dist = 1))
  }
  dev.off()
}
    



### --------------------------------
### plot celltype marker genes score 
tb = read.table('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/doc/Marker_Gene_MouseCorticalLayer_Soraia.csv',header=T,sep=',',as.is=T, row.names = 1)
tb = read.table('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/doc/Marker_Gene_Mouse_Soraia.csv',header=T,sep=',',as.is=T, row.names = 1)
tb = read.table('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/doc/Marker_Gene_Human_Soraia.csv',header=T,sep=',',as.is=T, row.names = 1)

af <- c('Marker_Gene_MouseCorticalLayer_Soraia.csv',
        'Marker_Gene_Mouse_Soraia.csv',
        'Marker_Gene_Human_Soraia.csv')
for (f in af){
  print(f)
  tb = read.table(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/doc/', f),header=T,sep=',',as.is=T, row.names = 1)
  for (i in seq(1,nrow(tb))){
    print(i)
    gvec = intersect(tb[i,], rownames(data))
    if (length(gvec) > 0){
      tmp <- t(sapply(gvec, function(g){
        tapply(data[g, ], list(clu), mean)
      }))
      clumean <- colMeans(tmp)
      V(mst)$color <- brewer.pal(9, 'Blues')[5]
      V(mst)$size <- clumean * 5
      png(paste0(pdir,paste0('MST_', sub('_.*','',sub('Marker_Gene_', '', f)), '_', rownames(tb)[i],'_markerScore.png')), width = 1800, height = 1800, res = 200)
      set.seed(12345)
      plot(mst, vertex.label = ct.clu, vertex.label.color="black", edge.color = 'black', edge.width=2, main = rownames(tb)[i])
      dev.off()
    }
  }
}
  

### ----------------------------
### plot marker genes on MSTtree
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





