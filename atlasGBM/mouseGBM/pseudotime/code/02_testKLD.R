library(here)
setwd(here())
library(Seurat)
library(LaplacesDemon)
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

pt.neuron <- seq(1, length(ord[[5]]))
names(pt.neuron) <- ord[[5]]
pt.neuron <- pt.neuron[grepl('GBM', names(pt.neuron))]
pt.oligo <- seq(1, length(ord[[16]]))
names(pt.oligo) <- ord[[16]]
pt.oligo <- pt.oligo[grepl('GBM', names(pt.oligo))]
pt.astrocyte <- seq(1, length(ord[[18]]))
names(pt.astrocyte) <- ord[[18]]
pt.astrocyte <- pt.astrocyte[grepl('GBM', names(pt.astrocyte))]

pairwiseKLD <- function(group1, group2){
  min = range(pt)[1]
  max = range(pt)[2]
  as.vector(sapply(group1, function(p1){
    sapply(group2, function(p2){
      a = density(pt[ap==p1], from = min, to = max)
      b = density(pt[ap==p2], from = min, to = max)
      KLD(a$y,b$y)$mean.sum.KLD
    })
  }))
}

groupKLD <- function(group1, group2){
  v1 = pairwiseKLD(group1, group2)
  v2 = as.vector(sapply(1:1000, function(i){
    if (!i%%10) print(i)
    set.seed(i)
    ap = c(group1, group2)
    gp1 = sample(ap, length(group1))
    gp2 = sample(ap, length(group2))
    mean(pairwiseKLD(gp1, gp2))  ## use mean only in permuation test
  }))
  ## KLD
  # KLD(density(v1, from = min(v1,v2), to = max(v1,v2))$y, density(v2, from = min(v1,v2), to = max(v1,v2))$y)$mean.sum.KLD
  ## t test
  # t.test(v1,v2)$p.value
  ## permutation test
  mean(v2 > mean(v1))
}

featureKLD <- function(feature){
  aval = unique(as.character(meta[,feature]))
  m = sapply(aval, function(val1){
    sapply(aval, function(val2){
      group1 <- rownames(meta[meta[, feature] == val1, , drop = FALSE])
      group2 <- rownames(meta[meta[, feature] == val2, , drop = FALSE])
      groupKLD(group1,group2)      
    })
  })
  m
}

for (traj in c('neuron', 'oligo','astrocyte')){
  if (traj == 'neuron')  pt <- pt.neuron
  if (traj == 'oligo') pt <- pt.oligo
  if (traj == 'astrocyte') pt <- pt.astrocyte
  ap <- sub('_.*', '', names(pt))
  loc <- unique(meta$Location)
  pdir <- paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/mouseGBM/pseudotime/plot/pt_', traj, '/')
  dir.create(pdir, showWarnings = F, recursive = T)
  rdir <- paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/mouseGBM/pseudotime/res/pt_', traj, '/')
  dir.create(rdir, showWarnings = F, recursive = T)
  
  min = min(pt)
  max = max(pt)
  up <- unique(ap)
  dmat <- matrix(NA,nrow=length(up),ncol=length(up),dimnames=list(up,up))
  for(p1 in 1:(length(up))){
    for (p2 in (p1):length(up)){
      a = density(pt[ap==up[p1]], from = min, to = max)
      b = density(pt[ap==up[p2]], from = min, to = max)
      dmat[up[p2],up[p1]] <- dmat[up[p1],up[p2]] <- KLD(a$y,b$y)$mean.sum.KLD
    }
  }
  write.csv(dmat, paste0(rdir, 'KLD_samples_pairwise.csv'))
  saveRDS(dmat, paste0(rdir, 'KLD_samples_pairwise.rds'))
  
  
  library(gplots)
  library(ggplot2)
  ## plot each sample pseudotime distribution
  feature = 'Location'
  pd <- data.frame(pseudotime = seq(1, length(pt)), 
                   sample = sub('_.*', '', names(pt)),
                   stringsAsFactors = F)
  pd <- pd[pd$sample %in% meta$study, ]
  pd[,feature] <- meta[match(pd$sample, meta$study), feature]
  pdf(paste0(pdir,feature,'_time_distribution.pdf'), width = 6, height = 7)
  ggplot(pd,aes_string(x=pd$pseudotime,fill=feature),alpha=0.5) + 
    geom_histogram() + 
    theme_classic() + 
    facet_wrap(~sample,scales = 'free_y',ncol=1,strip.position="right") + 
    theme(strip.text.y = element_text(angle = 0)) + 
          theme(axis.text=element_text(size=5))+
          scale_fill_brewer(palette='Set2')
  dev.off()
  
  ## plot each sample feature's significance heatmap (sig or insig)
  allfeat <- setdiff(colnames(meta), c('study', 'age', 'Pathology', 'Tumor.Grade', 'grade', 'Mutation.Rate'))
  for (feature in allfeat){
    print(feature)
    x=featureKLD(feature)
    write.csv(x, paste0(rdir, feature, '_significance.csv'))
    x[x<0.05] <- 0
    x[x>0.05] <- 1
    if (sum(x) != nrow(x)*ncol(x)){
      pdf(paste0(pdir,feature,'_significance.pdf'),width=6,height=6)
      par(mar=c(7,4,4,2)+0.1)
      print(heatmap.2(x, margins=c(12,8)))
      dev.off()  
    }
  }
}

  
