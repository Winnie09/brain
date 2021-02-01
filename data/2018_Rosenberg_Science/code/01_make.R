library(R.matlab)
d <- readMat('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Rosenberg_Science/raw/GSM3017261_150000_CNS_nuclei.mat')
mat <- matrix(0,nrow=nrow(d[[1]]),ncol=ncol(d[[1]]))
mat[1:50000,] <- as.integer(as.matrix(d[[1]][1:50000,]))
mat[50001:100000,] <- as.integer(as.matrix(d[[1]][50001:100000,]))
mat[100001:nrow(mat),] <- as.integer(as.matrix(d[[1]][100001:nrow(mat),]))
mat <- t(mat)
row.names(mat) <- d$genes
cn <- paste0('2018_Rosenberg:cell_',1:ncol(mat))
colnames(mat) <- cn
time <- d$sample.type
time[time %in% c('p11_brain','p11_spine')] <- 'postnatal day 11'
time[time %in% c('p2_brain ','p2_spine ')] <- 'postnatal day 2'
loc <- d$sample.type
loc[loc %in% c('p11_brain','p2_brain ')] <- 'brain'
loc[loc %in% c('p11_spine','p2_spine ')] <- 'spine'
ct <- d$cluster.assignment
ct <- sub(' .*','',ct)
ct <- as.numeric(ct)
ctv <- rep(NA,length(ct))
ctv[ct %in% 1:54] <- 'neuron'
ctv[ct %in% 55:60] <- 'oligodendrocyte'
ctv[ct %in% 61] <- 'oligodendroctye precursor cells'
ctv[ct %in% 62] <- 'perivivascular macrophage'
ctv[ct %in% 63] <- 'microglia'
ctv[ct %in% 64] <- 'endothelia cells'
ctv[ct %in% 65] <- 'smooth muscle cells'
ctv[ct %in% 66:67] <- 'vascular and leptomeningeal cells'
ctv[ct %in% 68:71] <- 'astrocytes'
ctv[ct %in% 72] <- 'ependyma'
ctv[ct %in% 73] <- 'schwann cells'
meta <- data.frame(cell=cn,celltype=ctv,time=time,location=loc,stringsAsFactors = F)
row.names(mat) <- gsub(' ','',row.names(mat))
mat <- mat[!duplicated(row.names(mat)),]
saveRDS(mat,file='/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Rosenberg_Science/proc/count.rds')
saveRDS(meta,file='/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Rosenberg_Science/proc/meta.rds')
