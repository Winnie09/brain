library(here)
setwd(here())
traj <- as.character(commandArgs(trailingOnly = T)[[1]])
i <- as.numeric(commandArgs(trailingOnly = T)[[2]])
j <- as.numeric(commandArgs(trailingOnly = T)[[3]])
print(traj)
print(i)
print(j)
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

if (traj == 'neuron') myord <- ord[[5]]
if (traj == 'oligo') myord <- ord[[16]]
if (traj == 'astrocyte') myord <- ord[[18]]

pt <- seq(1, length(myord))
names(pt) <- myord
ap <- sub('_.*', '', names(pt))
loc <- unique(meta$Location)
rdir <- paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/mouseGBM/pseudotime/res/testDensity/pt_', traj, '/', loc[i],'_',loc[j])
dir.create(rdir, showWarnings = F, recursive = T)

mat <- sapply(rownames(meta), function(p){
  v <- rep(0, length(myord))
  names(v) <- myord
  v[ap == p] <- 1
  v
}) ## cell by sample
mat.v <- as.vector(mat)  ## vectorize by columns
names(mat.v) <- paste0(rep(colnames(mat), each = nrow(mat)), ';', rownames(mat))
mat.m <- rbind(mat.v, mat.v+1)
rownames(mat.m) <- c('cellprop','cellpropFake')
cellanno <- data.frame(Cell = names(mat.v), Sample = sub(';.*', '', names(mat.v)), stringsAsFactors = F)

pseudotime = rep(pt, ncol(mat))
names(pseudotime) <- colnames(mat.m)

source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R')

group1 <- rownames(meta[meta$Location == loc[i], ])
group2 <- rownames(meta[meta$Location == loc[j], ])
ap.tmp <- c(group1, group2)
mat.tmp <- mat.m[, cellanno[,2] %in% ap.tmp]

cellanno.tmp <- cellanno[cellanno[,2] %in% ap.tmp, ]
design = data.frame(intercept = 1, location = ifelse(ap.tmp %in% group1, 0, 1), stringsAsFactors = FALSE)
rownames(design) <- ap.tmp
design <- as.matrix(design)
pt.tmp <- pseudotime[colnames(mat.tmp)]
res <- testpt(expr = mat.tmp, cellanno = cellanno.tmp, pseudotime = pt.tmp, design = design, test.type = 'Variable', ncores = 4)
saveRDS(res, paste0(rdir, '/testpt_density_res.rds'))

# res <- testpt(expr = mat.tmp, cellanno = cellanno.tmp, pseudotime = pt.tmp, design = design, type = 'Variable', ncores = 2, permuiter = 3)
# expr = mat.tmp; cellanno = cellanno.tmp; pseudotime = pt.tmp; design = design; type = 'Variable'; ncores = 1; permuiter = 3
# EMmaxiter=100; EMitercutoff=1; verbose=F; ncores=detectCores(); test.pattern = 'overall'; test.position = 'all'; fit.resolution = 1000; return.all.data = TRUE; demean = FALSE


