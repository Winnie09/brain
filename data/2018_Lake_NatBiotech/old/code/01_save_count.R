library(data.table)
d=fread('/scratch/users/whou10@jhu.edu/Wenpin/brain/data/2018_Lake_NatBiotech/raw/GSE97930_CerebellarHem_snDrop-seq_UMI_Count_Matrix_08-01-2017.txt',data.table=F)
rownames(d) = d[,1]
d = d[,-1]
d = as.matrix(d)
saveRDS(d, "/scratch/users/whou10@jhu.edu/Wenpin/brain/data/2018_Lake_NatBiotech/processed/count/CerebellarHem.rds")


d=fread('/scratch/users/whou10@jhu.edu/Wenpin/brain/data/2018_Lake_NatBiotech/raw/GSE97930_FrontalCortex_snDrop-seq_UMI_Count_Matrix_08-01-2017.txt',data.table=F)
rownames(d) = d[,1]
d = d[,-1]
d = as.matrix(d)
saveRDS(d, "/scratch/users/whou10@jhu.edu/Wenpin/brain/data/2018_Lake_NatBiotech/processed/count/FrontalCorte.rds")


d=fread('/scratch/users/whou10@jhu.edu/Wenpin/brain/data/2018_Lake_NatBiotech/raw/GSE97930_VisualCortex_snDrop-seq_UMI_Count_Matrix_08-01-2017.txt',data.table=F)
rownames(d) = d[,1]
d = d[,-1]
d = as.matrix(d)
saveRDS(d, "/scratch/users/whou10@jhu.edu/Wenpin/brain/data/2018_Lake_NatBiotech/processed/count/VisualCortex.rds")
