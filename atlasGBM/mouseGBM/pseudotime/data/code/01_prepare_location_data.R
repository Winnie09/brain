library(here)
setwd(here())
traj <- as.character(commandArgs(trailingOnly = T)[[1]])

pdir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/mouseGBM/pseudotime/plot/'

library(Seurat)
res <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/mouseGBM/integrate/harmony/36nonNormal_seuratGene/res/atlasGBM_harmony_with_ct.rds')
data = res@assays$RNA@data
harmony <- res@reductions$harmony@cell.embeddings
meta <- res@meta.data
meta <- meta[grepl('GBM', meta$study), 6:31]
meta <- meta[!duplicated(meta$study), ]
rownames(meta) <- meta$study
meta <- meta[,-c(3:10)]
colnames(meta)
meta <- meta[meta$Tumor.Grade == 'IV', ]

ord <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/mouseGBM/pseudotime/res/order_all.rds')

# > names(ord)
#  [1] "backbone 4,16,13,21,1,11,3" "branch: 4,5,19,25,23,9"    
#  [3] "branch: 4,16,35,14"         "branch: 4,16,35,15"        
#  [5] "branch: 4,16,13,21,1,12,18" "branch: 4,16,13,21,17,20"  
#  [7] "branch: 4,5,2,24"           "branch: 4,5,26"            
#  [9] "branch: 4,5,2,27"           "branch: 4,16,28"           
# [11] "branch: 4,5,29"             "branch: 4,16,13,21,1,30"   
# [13] "branch: 4,5,31"             "branch: 4,16,13,21,1,11,32"
# [15] "branch: 4,5,33"             "branch: 4,5,19,22,7,8,34"  
# [17] "branch: 4,5,6,36"           "branch: 4,10,37"   

path.neuron <- ord[[5]][grepl('GBM', ord[[5]])]
path.oligo <- ord[[16]][grepl('GBM', ord[[16]])]
path.astrocyte <- ord[[18]][grepl('GBM', ord[[18]])]
loc <- unique(meta$Location)

if (traj == 'pt_neuron') {
  ptorder.all <- path.neuron
  pseudotime.all <- seq(1, length(ord[[5]]))
  names(pseudotime.all) <- ord[[5]]
}
if (traj == 'pt_oligo') {
  ptorder.all <- path.oligo
  pseudotime.all <- seq(1, length(ord[[16]]))
  names(pseudotime.all) <- ord[[16]]
}
if (traj == 'pt_astrocyte') {
  ptorder.all <- path.astrocyte
  pseudotime.all <- seq(1, length(ord[[18]]))
  names(pseudotime.all) <- ord[[18]]
}

  
for (i in 1:3){
  print(i)
  for (j in (i+1):length(loc)){
    print(j)
    rdir <- paste0('atlasGBM/mouseGBM/pseudotime/data/testvar/',traj,'/', loc[i], '_', loc[j])
    dir.create(rdir, recursive = T, showWarnings = F)
    print(rdir)
    ## specify path, and samples of two locations
    ap = sub('_.*', '', ptorder.all)
    ptorder <- ptorder.all[ap %in% meta$study[meta$Location %in% c(loc[i], loc[j])]]
    ## prepare data
    expr <- as.matrix(data[, ptorder])
    cellanno = data.frame(Cell = colnames(expr), Sample = sub('_.*', '', colnames(expr)), stringsAsFactors = FALSE)
    pseudotime <- pseudotime.all[names(pseudotime.all) %in% colnames(expr)]
    design = data.frame(intercept = 1, location = ifelse(meta[as.character(unique(cellanno[,2])), 'Location'] == loc[i], 0, 1), stringsAsFactors = F)
    rownames(design) <- unique(cellanno[,2])
    design <- as.matrix(design)
    
    saveRDS(expr, paste0(rdir, '/expr.rds'))
    saveRDS(pseudotime, paste0(rdir, '/pseudotime.rds'))
    saveRDS(cellanno, paste0(rdir, '/cellanno.rds'))
    saveRDS(design, paste0(rdir, '/design.rds'))
  }
}



