library(here)
setwd(here())
pdir = '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/nmf/plot/pathwayscore_on_umap/'
pws = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/nmf/res/pathwayscore.rds')

umap = readRDS('atlasGBM/GBMonly/integrate/harmony/36nonNormalSeuratGene2000/res/umap_embeddings.rds')
umap = umap[rownames(pws), ]

pd = data.frame(u1 = umap[,1], u2 = umap[,2], pws)
library(ggplot2)
library(scattermore)
library(gridExtra)
library(RColorBrewer)

p1 <- ggplot(data = pd, aes(x=u1, y = u2, color = X1)) +
  geom_scattermore() +
  theme_classic()+
  scale_color_gradientn(colors = colorRampPalette(rev(brewer.pal(11,'RdYlBu')))(nrow(pd)))
p2 <- ggplot(data = pd, aes(x=u1, y = u2, color = X2)) +
  geom_scattermore() +
  theme_classic()+
  scale_color_gradientn(colors = colorRampPalette(rev(brewer.pal(11,'RdYlBu')))(nrow(pd)))
p3 <- ggplot(data = pd, aes(x=u1, y = u2, color = X3)) +
  geom_scattermore() +
  theme_classic()+
  scale_color_gradientn(colors = colorRampPalette(rev(brewer.pal(11,'RdYlBu')))(nrow(pd)))
p4 <- ggplot(data = pd, aes(x=u1, y = u2, color = X4)) +
  geom_scattermore() +
  theme_classic()+
  scale_color_gradientn(colors = colorRampPalette(rev(brewer.pal(11,'RdYlBu')))(nrow(pd)))

pdf(paste(pdir, 'GBM_integrate_umap.pdf'), width = 6, height = 4.5)
grid.arrange(p1,p2,p3,p4, ncol  = 2)
dev.off()


######3
umap = readRDS('atlasGBM/mouseGBM/integrate/harmony/36nonNormal_seuratGene/res/umap_embeddings.rds')
umap = umap[rownames(pws), ]
pd = data.frame(u1 = umap[,1], u2 = umap[,2], pws)
library(ggplot2)
library(scattermore)
library(gridExtra)
library(RColorBrewer)
p1 <- ggplot(data = pd, aes(x=u1, y = u2, color = X1)) +
  geom_scattermore() +
  theme_classic()+
  scale_color_gradientn(colors = colorRampPalette(rev(brewer.pal(11,'RdYlBu')))(nrow(pd))) 
p2 <- ggplot(data = pd, aes(x=u1, y = u2, color = X2)) +
  geom_scattermore() +
  theme_classic()+
  scale_color_gradientn(colors = colorRampPalette(rev(brewer.pal(11,'RdYlBu')))(nrow(pd)))
p3 <- ggplot(data = pd, aes(x=u1, y = u2, color = X3)) +
  geom_scattermore() +
  theme_classic()+
  scale_color_gradientn(colors = colorRampPalette(rev(brewer.pal(11,'RdYlBu')))(nrow(pd)))
p4 <- ggplot(data = pd, aes(x=u1, y = u2, color = X4)) +
  geom_scattermore() +
  theme_classic()+
  scale_color_gradientn(colors = colorRampPalette(rev(brewer.pal(11,'RdYlBu')))(nrow(pd)))

pdf(paste(pdir, 'mouseGBM_integrate_umap.pdf'), width = 6, height = 4.5)
grid.arrange(p1,p2,p3,p4, ncol  = 2)
dev.off()


#####
umap = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/humanGBM/integrate/harmony/36nonNormal_allen10xOnly_seuratGene/res/umap_embeddings.rds')
umap = umap[rownames(pws), ]
pd = data.frame(u1 = umap[,1], u2 = umap[,2], pws)
library(ggplot2)
library(scattermore)
library(gridExtra)
library(RColorBrewer)
p1 <- ggplot(data = pd, aes(x=u1, y = u2, color = X1)) +
  geom_scattermore() +
  theme_classic()+
  scale_color_gradientn(colors = colorRampPalette(rev(brewer.pal(11,'RdYlBu')))(nrow(pd))) 
p2 <- ggplot(data = pd, aes(x=u1, y = u2, color = X2)) +
  geom_scattermore() +
  theme_classic()+
  scale_color_gradientn(colors = colorRampPalette(rev(brewer.pal(11,'RdYlBu')))(nrow(pd)))
p3 <- ggplot(data = pd, aes(x=u1, y = u2, color = X3)) +
  geom_scattermore() +
  theme_classic()+
  scale_color_gradientn(colors = colorRampPalette(rev(brewer.pal(11,'RdYlBu')))(nrow(pd)))
p4 <- ggplot(data = pd, aes(x=u1, y = u2, color = X4)) +
  geom_scattermore() +
  theme_classic()+
  scale_color_gradientn(colors = colorRampPalette(rev(brewer.pal(11,'RdYlBu')))(nrow(pd)))

pdf(paste(pdir, 'humanGBM_integrate_umap.pdf'), width = 6, height = 4.5)
grid.arrange(p1,p2,p3,p4, ncol  = 2)
dev.off()


#########
allp = sub('_.*', '', rownames(pws))
pws.pb = t(sapply(unique(allp), function(p){
  tmp = pws[allp == p, , drop = F]
  colMeans(tmp)
}))
pr = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/pbpt/res/pca_expr.rds')
str(pr)

pr = pr[rownames(pws.pb),]
pd = data.frame(pc1 = pr[,1], pc2 = pr[,2], pws.pb)
library(ggplot2)
library(scattermore)
library(gridExtra)
library(RColorBrewer)
p1 <- ggplot(data = pd, aes(x=pc1, y = pc2, color = X1)) +
  geom_point() +
  theme_classic()+
  scale_color_gradientn(colors = colorRampPalette(rev(brewer.pal(11,'RdYlBu')))(nrow(pd))) 
p2 <- ggplot(data = pd, aes(x=pc1, y = pc2, color = X2)) +
  geom_point() +
  theme_classic()+
  scale_color_gradientn(colors = colorRampPalette(rev(brewer.pal(11,'RdYlBu')))(nrow(pd)))
p3 <- ggplot(data = pd, aes(x=pc1, y = pc2, color = X3)) +
  geom_point() +
  theme_classic()+
  scale_color_gradientn(colors = colorRampPalette(rev(brewer.pal(11,'RdYlBu')))(nrow(pd)))
p4 <- ggplot(data = pd, aes(x=pc1, y = pc2, color = X4)) +
  geom_point() +
  theme_classic()+
  scale_color_gradientn(colors = colorRampPalette(rev(brewer.pal(11,'RdYlBu')))(nrow(pd)))

pdf(paste(pdir, 'pb_pca_expr.pdf'), width = 6, height = 4.5)
grid.arrange(p1,p2,p3,p4, ncol  = 2)
dev.off()


####
pr = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/pbpt/res/pca_prop.rds')
str(pr)
allp = sub('_.*', '', rownames(pws))
pr = pr[rownames(pws.pb),]
pd = data.frame(pc1 = pr[,1], pc2 = pr[,2], pws.pb)
library(ggplot2)
library(scattermore)
library(gridExtra)
library(RColorBrewer)
p1 <- ggplot(data = pd, aes(x=pc1, y = pc2, color = X1)) +
  geom_point() +
  theme_classic()+
  scale_color_gradientn(colors = colorRampPalette(rev(brewer.pal(11,'RdYlBu')))(nrow(pd))) 
p2 <- ggplot(data = pd, aes(x=pc1, y = pc2, color = X2)) +
  geom_point() +
  theme_classic()+
  scale_color_gradientn(colors = colorRampPalette(rev(brewer.pal(11,'RdYlBu')))(nrow(pd)))
p3 <- ggplot(data = pd, aes(x=pc1, y = pc2, color = X3)) +
  geom_point() +
  theme_classic()+
  scale_color_gradientn(colors = colorRampPalette(rev(brewer.pal(11,'RdYlBu')))(nrow(pd)))
p4 <- ggplot(data = pd, aes(x=pc1, y = pc2, color = X4)) +
  geom_point() +
  theme_classic()+
  scale_color_gradientn(colors = colorRampPalette(rev(brewer.pal(11,'RdYlBu')))(nrow(pd)))

pdf(paste(pdir, 'pb_pca_prop.pdf'), width = 6, height = 4.5)
grid.arrange(p1,p2,p3,p4, ncol  = 2)
dev.off()

