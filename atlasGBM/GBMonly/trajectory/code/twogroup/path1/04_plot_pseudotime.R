ord = readRDS(file='/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/nmf/res/pseudotime_3_4.rds')
pt <- 1:length(ord)
names(pt) <- ord
facs <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/nmf/res/pathwayscore.rds')
facs = facs[names(pt), ]
meta <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/integrate/harmony/36nonNormalSeuratGene2000/res/meta.data.rds')
meta <- meta[,c('orig.ident','Location')]
meta <- unique(meta)
meta <- meta[meta[,1] %in% unique(sub('_.*','',ord)),]
uniloc = unique(meta[,2])
allp = sub('_.*', '', names(pt))

pdir = '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/trajectory/plot/twogroup/pseudotime/'
selcell = lapply(uniloc, function(i){
  names(pt)[allp %in% meta[meta[,2] %in% i,1]]
})
names(selcell) = uniloc

library(ggplot2)
library(RColorBrewer)
for (i in names(selcell)){
  print(i)
  pdf(paste0(pdir, 'pseudotime_3_4_', i,'.pdf'), width = 3.6, height = 2.8)
  mycolor = c(colorRampPalette(rev(brewer.pal(11, 'RdYlBu')))(length(selcell[[i]])), 'grey')
  df = data.frame(x = c(pd[!rownames(pd) %in% selcell[[i]],1], pd[selcell[[i]],1]), y = c(pd[!rownames(pd) %in% selcell[[i]], 2], pd[selcell[[i]],2]), time = c(rep(NA, sum(!rownames(pd) %in% selcell[[i]])), pt[selcell[[i]]]))
  set.seed(12345)
  df = df[sample(1:nrow(df), 1e5), ]
  print(ggplot(df) + 
          geom_point(aes(x = x, y = y, col = time), size = 0.1, stroke = 0) + 
          scale_color_gradientn(colors = mycolor) + 
          theme_classic() + xlab('factor score 3') + ylab('factor score 4') +
          ggtitle(i))
  dev.off()
}

pd = data.frame(x = facs[,3], y = facs[,4], pseudotime = pt, location = meta[match(allp, meta[,1]),2],stringsAsFactors = F)
pd[,4] = as.factor(pd[,4])
pdf(paste0(pdir, 'pseudotime_3_4.pdf'), width = 3.6, height = 2.3)
set.seed(12345)
df = pd[sample(1:nrow(pd), 1e5), ]
mycolor = c(colorRampPalette(rev(brewer.pal(11, 'RdYlBu')))(nrow(df)), 'grey')
print(ggplot(df) + 
        geom_point(aes(x = x, y = y, col = pseudotime), size = 0.1, stroke = 0) + 
        scale_color_gradientn(colors = mycolor) + 
        theme_classic() + xlab('factor score 3') + ylab('factor score 4'))
dev.off()

pdf(paste0(pdir, 'pseudotime_3_4_location.pdf'), width = 3.6, height = 2.3)
set.seed(12345)
df = pd[sample(1:nrow(pd), 1e5), ]
mycolor = c(colorRampPalette(rev(brewer.pal(11, 'RdYlBu')))(nrow(df)), 'grey')
ggplot(df) + 
        geom_point(aes(x = x, y = y, col = location), size = 0.1, stroke = 0) + 
        scale_color_brewer(palette = 'Dark2') + 
        theme_classic() + xlab('factor score 3') + ylab('factor score 4') +
  guides(color = guide_legend(override.aes = list(size = 3)))
dev.off()

