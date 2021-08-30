library(pheatmap)
library(Matrix)
path <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/nmf/res/overallgenepathway.rds')
ng <- nrow(path)
path <- as.vector(path)

expr <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/integrate/seurat/36nonNormalSeuratGene2000/res/norm.rds')
af <- list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cellcellcor/res',pattern = 'cellorder')

ord <- sapply(af, function(f){
  readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cellcellcor/res/',f))
})
ord.all = unlist(ord)

ord.gene = sapply(1:4, function(i){
  expr1 = expr[path[(1+30*(i-1)):(30*i)], ord.all]
  expr1 = expr1[, colSums(expr1)>0]
  pr = prcomp(expr1, scale. = T)$x
  cut = sort(cutree(hclust(dist(pr)), k = 5))
  names(cut)
})
ord.gene = as.vector(ord.gene)

for (f in af) {
  order <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cellcellcor/res/',f))
  pd <- expr[ord.gene,order]
  samp <- sub('_.*','',colnames(pd))
  for (s in unique(samp))
    pd[,samp==s] <- (pd[,samp==s]-rowMeans(pd[,samp==s]))
  png(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/nmf/plot/heatmap/',sub('cellorder_','',sub('.rds','',f)),'.png'),width=1000,height=1500,res=100)
  pheatmap(pd,cluster_rows = F,cluster_cols = F,show_colnames = F,color=colorRampPalette(c('blue','cyan','white','pink','red'))(20),gaps_col = cumsum(rle(sub('_.*','',colnames(pd)))$lengths),gaps_row = c(1:4)*30)
  dev.off()
}

