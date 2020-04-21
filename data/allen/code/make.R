meta = read.csv('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/allen/data/raw/sample_annotations.csv',as.is=T)
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/allen/data/raw')
library(scrattch.io)
options(stringsAsFactors = FALSE)
tome        <- "transcrip.tome"
exons       <- read_tome_dgCMatrix(tome,"data/t_exon")    # (or data/exon)
introns     <- read_tome_dgCMatrix(tome,"data/t_intron")  # (or data/intron)
sample_name <- read_tome_sample_names(tome)  
gene_name   <- read_tome_gene_names(tome)
# Note that the sample and gene names are NOT stored within the exon and intron matrices and must be read separately. 
m = exons + introns
mat <- matrix(0,nrow=50281,ncol=49494)
mat[,1:30000] <- as.integer(as.matrix(m[,1:30000]))
mat[,30001:ncol(m)] <- as.integer(as.matrix(m[,30001:ncol(m)]))
dimnames(mat) <- list(gene_name,sample_name)
id <- which(meta[,'class_label']!='Exclude')
mat <- mat[,id]
meta <- meta[id,]
cn <- paste0('allen:cell_',1:ncol(mat))
colnames(mat) <- cn
ct <- sub(' .*','',meta[,4])

conv <- c('Exc'='excitatory neurons','Inh'='inhibitory neurons','Gran'= 'granule cells','Purk1'= 'Purkinje neurons','Purk2'= 'Purkinje neurons', 'Endo'='endothelial cells', 'Astro'= 'astrocytes', 'Oligo'='oligodendrocytes', 'OPC'= 'oligodendrocytes precursor cells', 'Micro'='microglia','Peri'='pericytes','VLMC'='vascular and leptomeningeal cell')
layer <- sapply(meta[,'cluster_label'],function(i) strsplit(i,' ')[[1]][2],USE.NAMES = F)

meta <- data.frame(cell=cn,celltype=conv[ct],region=meta[,'region_label'],donor=meta[,'external_donor_name_label'],layer=layer,gender=meta[,'donor_sex_label'],stringsAsFactors = F)

saveRDS(meta,'/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/allen/data/proc/meta.rds')
saveRDS(mat, '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/allen/data/proc/count.rds')


