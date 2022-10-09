pdir <- '/home/whou10/scratch16/whou10/brain/predDNAm/predict/plot/'
rdir <- '/home/whou10/scratch16/whou10/brain/predDNAm/predict/res/'
source('/home/whou10/scratch16/whou10/resource/ggplot_theme.R')
theme_set(.new_theme)

d <- readRDS(paste0(rdir, 'cross_validation_correlation_columnSd.rds'))

library(ggplot2)
d2 <- do.call(rbind, d)


pd <- lapply(1:length(d), function(i){
  data.frame(d[[i]], test.data = names(d)[i])
})
pd <- do.call(rbind, pd)
br <- seq(min(pd[,2]), max(pd[,2]), length.out = 11)
br[2:(length(br)-1)] <- round(br[2:(length(br)-1)], 3)
c <- cut(pd[,2], breaks= br)
str(c)
length(c)

pd <- data.frame(pd, cut = c)
pd <- pd[complete.cases(pd), ]
dim(pd)
set.seed(123)
pd.tmp <- pd[sample(1:nrow(pd), 1e4, replace = FALSE), ]

pdf(paste0(pdir, 'pcc_by_sd.pdf'), width = 3.5, height = 1.8)
ggplot()+
  geom_violin(data = pd, aes(x = cut, y = correlation, fill = cut), scale = 'area', alpha = 0.2) +
  geom_jitter(data = pd.tmp, aes(x = cut, y = correlation, color = cut), size = 0.5, stroke = 0, alpha = 0.2, width = 0.2) + 
theme(axis.text.x = element_text(angle = 45, hjust = 0.5), legend.position = 'none')  +
  scale_fill_brewer(palette = 'RdYlBu', direction = -1)   + 
  scale_color_brewer(palette = 'RdYlBu', direction = -1)  + 
  xlab('Predicted value sd') + 
  ylab('Pearson correlation')
dev.off()


pdf(paste0(pdir, 'pcc_by_sd_stratify_testdata.pdf'), width = 8.5, height = 1.8)
ggplot()+
  geom_violin(data = pd, aes(x = cut, y = correlation, fill = test.data), scale = 'area', alpha = 0.2) +
  # geom_jitter(data = pd.tmp, aes(x = cut, y = correlation, color = cut), size = 0.5, stroke = 0, alpha = 0.2, width = 0.2) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5), legend.position = 'none')  +
  scale_fill_brewer(palette = 'RdYlBu', direction = -1)   + 
  scale_color_brewer(palette = 'RdYlBu', direction = -1)  + 
  xlab('Predicted value sd') + 
  ylab('Pearson correlation') +
  theme(legend.position = 'right')
dev.off()

