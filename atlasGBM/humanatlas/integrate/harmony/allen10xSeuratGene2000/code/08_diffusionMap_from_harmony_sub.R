library(Seurat)
library(destiny)
library(here)
here()
setwd(here())
rdir <- ddir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/humanatlas/integrate/harmony/allen10xSeuratGene2000/res/'
pdir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/humanatlas/integrate/harmony/allen10xSeuratGene2000/plot/diffusion_from_harmony_sub/'
dir.create(pdir, recursive = T)

# read in data, filtering
seu <- readRDS(paste0(ddir, 'humanAtlas_harmony.rds'))
ctAnno <- readRDS(paste0(ddir, 'celltype_annotation.rds'))
expr <- seu@assays$RNA@data ## librarysize-normalized  log2-transform
harmony <- seu@reductions$harmony@cell.embeddings
set.seed(12345)
harmony <- harmony[sample(1:nrow(harmony), 1e4), ]

# generate diffusion map
dm <- DiffusionMap(harmony) ## cell by harmany features
str(dm)
saveRDS(dm, paste0(rdir, 'diffusionMap.obj_from_pca_sub.rds'))
dmap <- dm@eigenvectors
rownames(dmap) <- rownames(harmony) ### !!!
str(dmap)
saveRDS(dmap, paste0(rdir, 'diffusionMap_from_pca_sub.rds'))

# plot diffusion map
pdf(paste0(pdir, '/dm_12dc.pdf'), width = 6, height = 5)
plot(dm,1:2,
pch = 20) 
dev.off()

pdf(paste0(pdir, '/dm_23dc.pdf'), width = 6, height = 5)
plot(dm,2:3,
pch = 20) 
dev.off()

pdf(paste0(pdir, '/dm_23dc.pdf'), width = 6, height = 5)
plot(dm,c(1,3),
pch = 20) 
dev.off()

pdf(paste0(pdir, '/dm_3dc.pdf'), width = 6, height = 6)
plot(dm,1:3,
pch = 20) 
dev.off()

# clustering
library(RColorBrewer)
library(scattermore)
library(ggplot2)
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/myfunc/01_function.R')
number.cluster <- NA
set.seed(12345)
clu <- mykmeans(t(dmap[,c(1,3)]), maxclunum = 50, number.cluster = number.cluster, ncores = 1, row.filter = FALSE)$cluster
table(clu)
saveRDS(clu, paste0(rdir, '/cluster.rds'))

pdf(paste0(pdir, '/cluster.pdf'), width = (0.8*max(clu))*2, height = (0.5*max(clu)))
p1 <- ggplot(data = data.frame(x = dmap[,1], y = dmap[,2], clu = as.factor(clu[rownames(dmap)])), 
             aes(x = x, y = y, color = clu)) +
  geom_scattermore()+
  scale_color_manual(values = colorRampPalette(brewer.pal(8,'Set1'))(max(clu)))+
  theme_classic() + xlab('DM1') + ylab('DM2')
p2 <- ggplot(data = data.frame(x = dmap[,1], y = dmap[,3], clu = as.factor(clu[rownames(dmap)])), 
             aes(x = x, y = y, color = clu)) +
  geom_scattermore()+
  scale_color_manual(values = colorRampPalette(brewer.pal(8,'Set1'))(max(clu)))+
  theme_classic() + xlab('DM1') + ylab('DM3')
gridExtra::grid.arrange(p1,p2, nrow=1)
dev.off()


pdf(paste0(pdir, '/dm_celltypeAnnotation.pdf'), width = 10, height = 4)
p1 <- ggplot(data = data.frame(x = dmap[,1], y = dmap[,2], ct = as.factor(ctAnno[rownames(dmap)])), aes(x = x, y = y, color = ct)) +
  geom_scattermore()+
  scale_color_manual(values = colorRampPalette(brewer.pal(8,'Set1'))(length(unique(ctAnno))))+
  theme_classic() + xlab('DM1') + ylab('DM2')
p2 <- ggplot(data = data.frame(x = dmap[,1], y = dmap[,3], ct = as.factor(ctAnno[rownames(dmap)])), 
             aes(x = x, y = y, color = ct)) +
  geom_scattermore()+
    scale_color_manual(values = colorRampPalette(brewer.pal(8,'Set1'))(length(unique(ctAnno))))+
  theme_classic() + xlab('DM1') + ylab('DM3')
gridExtra::grid.arrange(p1,p2, nrow=1)
dev.off()


