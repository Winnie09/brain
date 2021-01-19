library(destiny)
# study <- '2018_Li_Science'
study <- as.character(commandArgs(trailingOnly = TRUE)[[1]])
print(study)
ddir <- paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/', study, '/pca/')
pdir <- rdir <- paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/', study, '/dm/')
dir.create(rdir, showWarnings = FALSE, recursive = TRUE)
# read in data, filtering
pr <- readRDS(paste0(ddir, 'pr.rds'))

# generate diffusion map
dm <- DiffusionMap(pr)
str(dm)
saveRDS(dm, paste0(rdir, 'diffusionMap.obj.rds'))
dmap <- dm@eigenvectors
rownames(dmap) <- rownames(pr)
str(dmap)
saveRDS(dmap, paste0(rdir, 'diffusionMap.rds'))

# plot diffusion map
pdf(paste0(pdir, '/dm_12dc.pdf'), width = 6, height = 5)
plot(dm,1:2,
pch = 20) 
dev.off()

pdf(paste0(pdir, '/dm_13dc.pdf'), width = 6, height = 5)
plot(dm,1:3,
pch = 20) 
dev.off()

pdf(paste0(pdir, '/dm_23dc.pdf'), width = 6, height = 5)
plot(dm,2:3,
pch = 20) 
dev.off()

