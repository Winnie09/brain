library(Seurat)
library(destiny)
library(here)
here()
setwd(here())
rdir <- ddir <- 'atlasGBM/subHumanGBM/integrate/harmony/rmCt/result/'
pdir <- 'atlasGBM/subHumanGBM/integrate/harmony/rmCt/plot/'

# read in data, filtering
seu <- readRDS(paste0(ddir, 'atlasGBM_harmony.rds'))
expr <- seu@assays$RNA@data ## librarysize-normalized  log2-transform

len <- ncol(expr)/1e4 + 1
exprlist <- lapply(1:len, function(i){
  as.matrix(expr[, seq(1e4*(i-1) + 1, min(1e4*i, ncol(expr)))])
})
expr <- do.call(cbind, exprlist)
str(expr)
expr <- expr[rowMeans(expr>0.1)>0.01, ] ## nothing changed
str(expr) ##

# generate diffusion map
dm <- DiffusionMap(t(expr), )
str(dm)
saveRDS(dm, paste0(rdir, 'diffusionMap.obj.rds'))
dmap <- dm@eigenvectors
rownames(dmap) <- colnames(expr)
str(dmap)
saveRDS(dmap, paste0(rdir, 'diffusionMap.rds'))

# plot diffusion map
pdf(paste0(plotdir, '/dm_12dc.pdf'), width = 6, height = 5)
plot(dm,1:2,
pch = 20) 
dev.off()

pdf(paste0(plotdir, '/dm_23dc.pdf'), width = 6, height = 5)
plot(dm,2:3,
pch = 20) 
dev.off()

pdf(paste0(plotdir, '/dm_23dc.pdf'), width = 6, height = 5)
plot(dm,c(1,3),
pch = 20) 
dev.off()

pdf(paste0(plotdir, '/dm_3dc.pdf'), width = 6, height = 6)
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

