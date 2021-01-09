setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/')     
rdir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/humanatlas/data/'
mat <- readRDS(paste0(rdir, 'combine_mat.rds'))

meta <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/humanGBM/data/meta_allcell.rds')
# rownames(meta) <- sub('.*;', '', meta$cell)
rownames(meta) <- meta$cell
a = meta[,c(-1,-2,-6,-7,-32)]  
meta.a <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/allen_human_10x/proc/meta.rds')

## change meta.a cell names and rownames
v <- meta.a$cell
v1 <- sub(':.*', '', v)
v2 <- sub('.*:', '', v)
v.new <- paste0('allen_human_10x;allen_human_10x:', v2)
meta.a$cell <- v.new
rownames(meta.a) <- v.new
## change meta cell names and rownames
v <- meta$cell
id <- grep('allen', v)
tmp <- v[id]
v1 <- sub(':.*', '', tmp)
v2 <- sub('.*:', '', tmp)
v[id] <- paste0('allen_human_fluidigm;allen_human_fluidigm:',v2)
meta$cell <- v
rownames(meta) <- v


## subset meta
# exist.cell <- sub('.*;', '', meta$cell)
exist.cell <- meta$cell
int = intersect(colnames(mat), exist.cell)
exist.meta <- meta[!grepl('GBM',meta$study), ]
# rownames(exist.meta) <- sub('.*;', '', exist.meta$cell)
rownames(exist.meta) <- exist.meta$cell
exist.meta$study <- sapply(exist.meta$study, function(i) ifelse(i == 'allen', 'allen_human_fluidigm', i))

## make new meta
identical(colnames(mat), c(rownames(exist.meta)))
meta.a$gender <- ifelse(meta.a$gender == 'F', 'female', 'male')
meta.a$study <- 'allen_human_10x'
meta.a$species <- 'human'
meta.a$cell.shortname <- sub('.*;', '', meta.a$cell)

tmp <- matrix(NA, nrow = nrow(meta.a), ncol = length(setdiff(colnames(exist.meta), colnames(meta.a))))
tmp <- as.data.frame(tmp, stringsAsFactors = FALSE)
colnames(tmp) <- setdiff(colnames(exist.meta), colnames(meta.a))
meta.a.new <- cbind(meta.a, tmp)

meta.new <- rbind(exist.meta, meta.a.new[, colnames(exist.meta)])

length(setdiff(colnames(mat), meta.new$cell.shortname))
length(setdiff(colnames(mat), sub('.*;', '', meta.new$cell)))

meta.new2 <- meta.new[match(colnames(mat),  meta.new$cell), ]
identical(colnames(mat), meta.new2$cell)

meta.new2$celltype[meta.new2$celltype == 'astrocytes'] <- 'astrocyte'
meta.new2$celltype[meta.new2$celltype == 'endothelial cells'] <- 'endothelial'
meta.new2$celltype[meta.new2$celltype == 'endothelials'] <- 'endothelial'
meta.new2$celltype[meta.new2$celltype == 'oligodendrocytes precursor cells'] <- 'oligodendrocytes precursor'
meta.new2$celltype[meta.new2$celltype == 'oligodendrocytes'] <- 'oligodendrocyte'
meta.new2$celltype[meta.new2$celltype == 'vascular and leptomeningeal cell'] <- 'vascular and leptomeningeal'

saveRDS(meta.new2, paste0(rdir, 'meta_allcell.rds'))

