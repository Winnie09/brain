library(Seurat)
library(destiny)
library(here)
here()
setwd(here())
rdir <- ddir <- 'atlasGBM/subMouseGBM/integrate/harmony/rmCt_onlyGradeIV/result/'
pdir <- 'atlasGBM/subMouseGBM/integrate/harmony/rmCt_onlyGradeIV/plot/diffusion_from_pca_sub/'
dir.create(pdir, recursive = T)

# read in data, filtering
seu <- readRDS(paste0(ddir, 'atlasGBM_harmony.rds'))
expr <- seu@assays$RNA@data ## librarysize-normalized  log2-transform
harmony <- seu@reductions$harmony@cell.embeddings
set.seed(12345)
harmony <- harmony[sample(1:nrow(harmony), 1e4), ]

# generate diffusion map
dm <- DiffusionMap(harmony) ## cell by harmany features
str(dm)
saveRDS(dm, paste0(rdir, 'diffusionMap.obj_from_pca_sub.rds'))
dmap <- dm@eigenvectors
# rownames(dmap) <- colnames(expr)
rownames(dmap) <- rownames(harmony)
str(dmap)
saveRDS(dmap, paste0(rdir, 'diffusionMap_from_pca_sub.rds'))



# plot diffusion map
ctanno.atlas <- readRDS('atlasGBM/mouseatlas/integrate/harmony/seuratGene/res/celltype_annotation.rds')
ctanno.gbm <- readRDS('atlasGBM/GBMonly/integrate/harmony/36nonNormalSeuratGene2000/res/celltype_annotation.rds')
ctanno <- c(ctanno.atlas, ctanno.gbm)


pd <- reshape2::melt(dmap)
colnames(pd) <- c('cell', 'DC', 'expr')
pd <- cbind(pd, ctanno[match(as.character(pd[,1]), names(ctanno))])
pdf(paste0(plotdir, '/dm_12dc.pdf'), width = 6, height = 5)
plot(dmap,1:2,
pch = 20) 
dev.off()

pdf(paste0(plotdir, '/dm_23dc.pdf'), width = 6, height = 5)
plot(dmap,2:3,
pch = 20) 
dev.off()

pdf(paste0(plotdir, '/dm_23dc.pdf'), width = 6, height = 5)
plot(dmap,c(1,3),
pch = 20) 
dev.off()

pdf(paste0(plotdir, '/dm_3dc.pdf'), width = 6, height = 6)
plot(dmap,1:3,
pch = 20) 
dev.off()

# clustering
library(RColorBrewer)
library(scattermore)
library(ggplot2)
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/myfunc/01_function.R')
number.cluster <- NA
set.seed(12345)
clu <- mykmeans(dmap[,c(1,3)], maxclunum = 50, number.cluster = number.cluster)$cluster
table(clu)
saveRDS(clu, paste0(rdir, '/cluster.rds'))
pdf(paste0(plotdir, '/cluster.pdf'), width = (0.8*max(clu))*2, height = (0.5*max(clu)))
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






