library(data.table)
library(mixtools)
library(parallel)
af <- list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_myelymreference/',pattern = '^GBM')
dum <- mclapply(af,function(f) {
  print(f)
  #library(gplots)
  #image(d,col=bluered(101))
  a <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/cnv_myelymreference/',f,'/cutoff0.1/output/run.final.infercnv_obj'))
  chr <- a@gene_order
  d <- a@expr.data
  d <- d - 1
  chr <- chr[rownames(d),1]
  
  refd <- d[,a@tumor_subclusters$hc$all_references$labels]
  
  d <- d[,a@tumor_subclusters$hc$all_observations$labels]
  d <- d[,a@tumor_subclusters$hc$all_observations$order]
  
  bd <- matrix(0,nrow=nrow(d),ncol=ncol(d),dimnames = dimnames(d))
  
  k <- refd[refd < -0.01]
  negcut <- mean(k)-1*sd(k)
  
  k <- refd[refd > 0.01]
  poscut <- mean(k)+1*sd(k)  
  
  bd[d > poscut] <- 1
  bd[d < negcut] <- -1
  
  cl <- sapply(unique(chr),function(i) {
    id <- which(chr==i)
    winlen <- round(length(id)/3)
    posprop <- sapply(1:(length(id)-winlen),function(j) {
      colMeans(bd[id[j:(j+winlen)],]==1) > 0.8
    })
    negprop <- sapply(1:(length(id)-winlen),function(j) {
      colMeans(bd[id[j:(j+winlen)],]==-1) > 0.8
    })
    list(pos=names(which(rowSums(posprop) > 1)),neg=names(which(rowSums(negprop) > 1)))
  },simplify = F)
  names(cl) <- unique(chr)
  
  # pd <- matrix(0,nrow=length(cl),ncol=ncol(d),dimnames = list(names(cl),colnames(d)))
  # for (i in names(cl)) {
  #   pd[i,cl[[i]][['pos']]] <- 1
  #   pd[i,cl[[i]][['neg']]] <- -1
  # }
  # pdf(paste0('/home-4/zji4@jhu.edu/scratch/cnvsel/clu/',f,'.pdf'))  
  # image(pd)
  # dev.off()
  saveRDS(cl,file=paste0('/home-4/zji4@jhu.edu/scratch/cnvsel/cnv/',f,'.rds'))
},mc.cores=detectCores())



