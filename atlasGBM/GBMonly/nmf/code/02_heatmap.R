library(pheatmap)
path <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/nmf/res/overallgenepathway.rds')
ng <- nrow(path)
path <- as.vector(path)

expr <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/data/36nonNormal_combined_log2norm_mat.rds')

af <- list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cellcellcor/res',pattern = 'cellorder')
for (f in af) {
  order <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cellcellcor/res/',f))
  
  pd <- expr[path,order]
  samp <- sub('_.*','',colnames(pd))
  for (s in unique(samp))
    pd[,samp==s] <- (pd[,samp==s]-rowMeans(pd[,samp==s]))
  
  # order <- unlist(sapply(unique(samp),function(s) {
  #   colnames(pd)[samp==s][hclust(dist(t(pd[,samp==s])))$order]
  # }))
  # pd <- pd[,order]
  
  # pd[pd > quantile(pd,0.95)] <- quantile(pd,0.95)
  # pd[pd < quantile(pd,0.05)] <- quantile(pd,0.05)
  png(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/nmf/plot/heatmap/',sub('cellorder_','',sub('.rds','',f)),'.png'),width=1000,height=1500,res=100)
  pheatmap(pd,cluster_rows = F,cluster_cols = F,show_colnames = F,color=colorRampPalette(c('blue','cyan','white','pink','red'))(20),gaps_col = cumsum(rle(sub('_.*','',colnames(pd)))$lengths),gaps_row = c(1:4)*30)
  dev.off()  
}

