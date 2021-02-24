library(here)
setwd(here())

meta <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/mouseGBM/integrate/harmony/36nonNormal_seuratGene/res/meta.data.rds')
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


traj = 'neuron'
if (traj == 'neuron') myord <- ord[[5]]
if (traj == 'oligo') myord <- ord[[16]]
if (traj == 'astrocyte') myord <- ord[[18]]

pt <- seq(1, length(myord))
names(pt) <- myord
ap <- sub('_.*', '', names(pt))

loc <- unique(meta$Location)
pdir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/mouseGBM/pseudotime/plot/pt_neuron/'

rdir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/mouseGBM/pseudotime/res/pt_neuron/'

mat <- sapply(rownames(meta), function(p){
  v <- rep(0, length(myord))
  names(v) <- myord
  v[ap == p] <- 1
  v
}) ## cell by sample
mat.v <- as.vector(mat)  ## vectorize by columns
names(mat.v) <- paste0(rep(colnames(mat), each = nrow(mat)), ';', rownames(mat))
mat.m <- t(as.matrix(mat.v, nrow = 1))
rownames(mat.m) <- 'cellprop'

  
cellanno <- data.frame(Cell = names(mat.v), Sample = sub(';.*', '', names(mat.v)), stringsAsFactors = F)

pseudotime = rep(pt, ncol(mat))
names(pseudotime) <- colnames(mat.m)


design = data.frame(intercept = 1, location = meta$Location, stringsAsFactors = F)
rownames(design) <- rownames(meta)

source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R')

res <- testpt(expr = mat.m, cellanno = cellanno, pseudotime = pseudotime, design = design, type = 'Variable', ncores = 10)
saveRDS(res, 'testpt_density_res.rds')
