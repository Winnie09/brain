library(loomR)
lfile <- connect(filename = "/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Zeisel_Cell/raw/l5_all.loom", mode = "r+")
#lfile <- connect(filename = "/home-4/zji4@jhu.edu/scratch/brain/data/2018_Zeisel_Cell/raw/l5_all.loom", mode = "r+")
mat <- lfile[['matrix']][,]
gn <- lfile[["row_attrs/gene_names"]][]
cn <- lfile[["col_attrs/gene_names"]][]
rownames(mat) = gn
colnames(mat) = cn
row.names(mat) <- gsub(' ','',row.names(mat))
saveRDS(mat,'/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Zeisel_Cell/proc/mat.rds')

