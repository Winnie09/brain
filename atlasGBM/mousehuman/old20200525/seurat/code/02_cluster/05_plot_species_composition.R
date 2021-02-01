setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlas/seurat/result')
# key = as.character(commandArgs(trailingOnly = T)[[1]])
key = 'atlasOnly_homolog_sub'
# key = 'atlas_GBM_homolog_sub'
library(Seurat)
clu = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlas/seurat/result/atlasOnly_homolog_sub/k3_23cluster/cluster/cluster.rds')
meta = readRDS(paste0('./',key,'/meta.rds'))
df <- data.frame(clu = clu, cell = names(clu), species = as.character(meta[match(names(clu), meta$cell), 'species']))

library(ggplot2)
pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlas/seurat/plot/atlasOnly_homolog_sub/umap_cluster23/species_composition.pdf',width=6,height=6)
ggplot(data=df, aes(x=clu, fill=species)) +
  geom_density()+
  scale_fill_brewer(palette="Paired")+
  theme_minimal() +
  xlab('Cluster') + ylab('Species Composition (Num.cells)')
dev.off()



res1 <- do.call(rbind,lapply(seq(1,length(unique(df[,1]))), function(i){
  a <- nrow(df[df[,1]==i, ])
  b <- nrow(df[df[,1]==i & df[,3]=='human', ])
  data.frame(i, round(b/a,3),species='human',stringsAsFactors = F)
}))
res2 <- do.call(rbind,lapply(seq(1,length(unique(df[,1]))), function(i){
  a <- nrow(df[df[,1]==i, ])
  b <- nrow(df[df[,1]==i & df[,3]=='mouse', ])
  data.frame(i, round(b/a,3),species='mouse',stringsAsFactors = F)
}))
colnames(res1) <- colnames(res2) <- c('cluster','proportion','species')
res <- rbind(res1,res2)
res[,1] <- factor(res[,1])
pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlas/seurat/plot/atlasOnly_homolog_sub/umap_cluster23/species_composition_prop.pdf',width=8,height=4)
ggplot(res,aes(x=cluster,y=proportion,fill=species)) + geom_bar(position="dodge", stat="identity") + 
  scale_fill_brewer(palette="Paired")+
  theme_minimal() +
  xlab('Cluster') + ylab('Species Composition (Proportion)')
dev.off()


