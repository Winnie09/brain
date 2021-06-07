library(RColorBrewer)
pathwayscore = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/nmf/res/pathwayscore.rds')

pdir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/nmf/plot/pathwayscore/celllcellcor_order/'
dir.create(pdir)
# co <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cellcellcorsinglecenter/res/cellorder_Frontal.rds')
af <- list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cellcellcor/res/')
af <- af[grepl('cellorder_', af)]
for (f in af){
  print(f)
  co <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cellcellcor/res/', f))
  mat = pathwayscore[co, ]
  allp <- sub('_.*', '', rownames(mat))
  for (p in unique(allp)){
    png(paste0(pdir, p,'.png'), width = 300, height = 300, res = 100)
    heatmap(mat[allp == p, ], scale = 'none', Rowv = NA, Colv = NA, labRow = NA, labCol = NA,col=colorRampPalette(c('blue3','blue','lightblue','white','pink','red','red3'))(100))
    dev.off()
  }
}

####
pdir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/nmf/plot/pathwayscore/celllcellcorsinglecenter_order/'
dir.create(pdir)
af <- list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cellcellcorsinglecenter/res/')
af <- af[grepl('cellorder_', af)]
for (f in af){
  print(f)
  co <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cellcellcorsinglecenter/res/', f))
  mat = pathwayscore[co, ]
  allp <- sub('_.*', '', rownames(mat))
  for (p in unique(allp)){
    png(paste0(pdir, p,'.png'), width = 300, height = 300, res = 100)
    heatmap(mat[allp == p, ], scale = 'none', Rowv = NA, Colv = NA, labRow = NA, labCol = NA,col=colorRampPalette(c('blue3','blue','lightblue','white','pink','red','red3'))(100))
    dev.off()
  }
}

