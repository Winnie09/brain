library(pheatmap)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(ape)
pdir = '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/nmf/plot/pathwayscore/cellcellcor_order_cnvclu/'
dir.create(pdir, showWarnings = F)
## read cell cell correlation order for a sample
allf = list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cellcellcor/res/', pattern = 'cellorder_')
for (f in allf){
  print(f)
  co <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cellcellcor/res/', f))
  allp <- sub('_.*', '', co)
  for (p in unique(allp)){
    print(p)
    ## get cnv clusters and maglinancy clusters
    library(ape)
    d <- read.tree(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_allenreference/',p,'/cutoff0.1/output/infercnv.observations_dendrogram.txt'))
    clu <- cutree(as.hclust(d),10)
    cos = rev(co[allp == p])
    cos = intersect(cos, names(clu))
    
    myclu = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_allenreference/res/cnvchrclu/', p, '.rds'))
    mat = matrix(c(clu[cos], myclu[cos]), ncol = 2)
    dimnames(mat) = list(cos, c('cnvclu', 'cnvchrclu'))
    ## read pathwayscore
    pas = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/nmf/res/pathwayscore.rds')
    pas = pas[cos,]
    colnames(pas) = paste0('pathway',1:4)
    pas = apply(pas, 2, scale)
    ht1 = Heatmap(pas,cluster_rows = F,
                  cluster_columns = FALSE,
                  show_row_names = F,
                  show_column_names = T,
                  col = colorRampPalette(c('blue3','blue','lightblue','white','pink','red','red3'))(100))
    ht2 = Heatmap(mat, 
                  cluster_rows = F,
                  cluster_columns = FALSE,
                  show_row_names = F,
                  show_column_names = T,
                  col = colorRampPalette(brewer.pal(n = 8, name ="Set1"))(length(unique(mat[,1]))+5)[-c(2,4,9,13,14)])
    htlist <- ht1 + ht2 
    pdf(paste0(pdir, p, '.pdf'), width = 5.3, height = 3.3)
    draw(htlist, merge_legend = FALSE,annotation_legend_side = "right",heatmap_legend_side = "right")
    dev.off()
  }
}

