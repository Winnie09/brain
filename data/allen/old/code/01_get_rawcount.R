setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/allenBrain/data/raw/')
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
sum((colSums(mat)-colSums(m))^2)
saveRDS(mat, '/home-4/whou10@jhu.edu/scratch/Wenpin/allenBrain/data/proc/matrix/rawcount.rds')
rc <- colSums(mat)
gn <- colSums(mat > 0)
mito <- colSums(mat[grep('^MT-',row.names(mat),ignore.case=T),])/rc
savedir = '/home-4/whou10@jhu.edu/scratch/Wenpin/allenBrain/data/proc/qc/'
dir.create(savedir,showWarnings = F)
saveRDS(rc,file=paste0(savedir,'totalreadcount.rds'))
saveRDS(gn,file=paste0(savedir,'expressedgenenumber.rds'))
saveRDS(mito,file=paste0(savedir,'mitoproportion.rds'))

