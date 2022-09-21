# Cells were classified as cancer and normal in previous analysis in our GBM scRNA-seq data. Gene expression raw counts were aggregated across cancer or normal cells within each sample to create pseudobulks, followed by library size normalization and log2 transformation after adding a pseudocount of 1. Limma was used to identify the differential genes between cancer and normal pseudobulks pairing each sample. logFC (base 2) indicates the difference between cancer and normal. Positive logFC means the genes have higher expression in cancer cells versus normal cells. 

mat <- readRDS('/home/whou10/scratch16/whou10/brain/atlasGBM/GBMonly/data/36nonNormal_combined_count_mat.rds')
cancercell <- readRDS('/home/whou10/scratch16/whou10/brain/atlasGBM/GBMonly/cnv_myelymreference/res/cancercell_final.rds')
normalcell <- readRDS('/home/whou10/scratch16/whou10/brain/atlasGBM/GBMonly/cnv_myelymreference/res/normalcell_final.rds')

canmat <- mat[,names(cancercell)]
normat <- mat[,names(normalcell)]
canmat <- rowsum(t(canmat),sub('_.*','',colnames(canmat)))
normat <- rowsum(t(normat),sub('_.*','',colnames(normat)))

canmat <- t(log2(canmat/rowSums(canmat) * 1e7 + 1))
normat <- t(log2(normat/rowSums(normat) * 1e7 + 1))
difmat <- canmat-normat[,colnames(canmat)]

library(limma)
res <- topTable(eBayes(lmFit(difmat,matrix(1,nrow=ncol(difmat)))),n=nrow(difmat))
res <- res[order(-res[,1]),]

write.csv(res,file='/home/whou10/scratch16/whou10/brain/atlasGBM/GBMonly/cancer_vs_normal/diff.csv')

res2 = res[res[,1] > 0 & res[,5] < 0.05, ]
str(res2) # > 11,000
res3 = res[res[,1] > 1 & res[,5] < 0.05, ]
str(res3) # 3658

write.csv(res3, file = '/home/whou10/scratch16/whou10/brain/atlasGBM/GBMonly/cancer_vs_normal/diffgene_high_in_cancer_logFC1_FDR0.05.csv')


