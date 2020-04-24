m = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/allenBrain/data/proc/matrix/count.rds')
meta = read.csv('/home-4/whou10@jhu.edu/scratch/Wenpin/allenBrain/data/raw/sample_annotations.csv')
ct = meta[match(colnames(m),meta$sample_name), 'cell_type_alias_label']
ct = sub(' .*','',ct)
names(ct) = colnames(m)
saveRDS(ct,'/home-4/whou10@jhu.edu/scratch/Wenpin/allenBrain/data/proc/matrix/celltype.rds')



