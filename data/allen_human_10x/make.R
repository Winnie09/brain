meta = read.csv('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/allen_human_10x/raw/metadata.csv',as.is=T)
library(data.table)
d <- fread('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/allen_human_10x/raw/matrix.csv',data.table=F)
rownames(d) <- d[,1]
d <- as.matrix(d[,-1])
d <- t(d)

cn <- paste0('allen:cell_',1:ncol(d))
colnames(d) <- cn
ct <- sub(' .*','',meta[,'cluster_label'])

conv <- c('Exc'='excitatory neurons','Inh'='inhibitory neurons','Gran'= 'granule cells','Purk1'= 'Purkinje neurons','Purk2'= 'Purkinje neurons', 'Endo'='endothelial cells', 'Astro'= 'astrocytes', 'Oligo'='oligodendrocytes', 'OPC'= 'oligodendrocytes precursor cells', 'Micro'='microglia','Peri'='pericytes','VLMC'='vascular and leptomeningeal cell')
layer <- sapply(meta[,'cluster_label'],function(i) strsplit(i,' ')[[1]][2],USE.NAMES = F)

meta <- data.frame(cell=cn,celltype=conv[ct],region=meta[,'region_label'],donor=meta[,'external_donor_name_label'],layer=layer,gender=meta[,'donor_sex_label'],stringsAsFactors = F)

saveRDS(meta,'/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/allen_human_10x/proc/meta.rds')
saveRDS(d, '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/allen_human_10x/proc/count.rds')
