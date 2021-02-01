load('/scratch/users/whou10@jhu.edu/Wenpin/brain/data/2018_Li_Science/raw/Sestan.adultHumanNuclei.Psychencode.Rdata')
row.names(umi2) <- sub('.*\\|','',row.names(umi2))
umi2 <- umi2[!duplicated(row.names(umi2)),]
cn <- paste0('2018_Li_Science:adult:cell_',1:ncol(umi2))
ctv <- meta2$ctype
conv <- c('ExN'='excitatory neurons','InN'='inhibitory neurons','Gran'= 'granule cells','Purk1'= 'Purkinje neurons','Purk2'= 'Purkinje neurons', 'Endo'='endothelial cells', 'Astro'= 'astrocytes', 'Oligo'='oligodendrocytes', 'OPC'= 'oligodendrocytes precursor cells', 'Microglia'='microglia','Per'='pericytes','VSMC'='vascular smooth muscle cells')
m <- data.frame(cell=cn,celltype=conv[ctv],time='adult',species='human',stringsAsFactors = F)
colnames(umi2) <- cn
saveRDS(umi2,file='/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Li_Science/proc/expr/adult.rds')
saveRDS(m,file='/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Li_Science/proc/meta_adult.rds')

load('/scratch/users/whou10@jhu.edu/Wenpin/brain/data/2018_Li_Science/raw/Sestan.fetalHuman.Psychencode.Rdata')
umi2 <- as.matrix(count2)
row.names(umi2) <- sub('.*\\|','',row.names(umi2))
umi2 <- umi2[!duplicated(row.names(umi2)),]
cn <- paste0('2018_Li_Science:fetal:cell_',1:ncol(umi2))
ctv <- meta2$ctype
cctv <- ctv
cctv[grep('ExN',cctv)] <- 'excitatory neurons'
cctv[grep('InN',cctv)] <- 'inhibitory neurons'
cctv[grep('NasN',cctv)] <- 'nascent neurons'
cctv[grep('OPC',cctv)] <- 'oligodendrocytes precursor cells'
cctv[grep('Astro',cctv)] <- 'astrocytes'
cctv[grep('Endo',cctv)] <- 'endothelial cells'
cctv[grep('Microglia',cctv)] <- 'microglia'
cctv[grep('Oligo',cctv)] <- 'oligodendrocytes'
cctv[grep('Pericyte',cctv)] <- 'pericytes'
cctv[grep('NEPRGC',cctv)] <- 'neural epithelial progenitor/radial glial cells'
cctv[grep('IPC',cctv)] <- 'intermediate progenitor cells'

m1 <- data.frame(cell=cn,celltype=cctv,gender=meta2[,'Sex'],time=sub('PCW',' post-conception weeks',meta2[,'Age']),species='human',stringsAsFactors = F)
colnames(umi2) <- cn
saveRDS(umi2,file='/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Li_Science/proc/expr/fetal.rds')
saveRDS(m1,file='/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Li_Science/proc/meta_fetal.rds')

m$gender <- 'unknown'
meta = rbind(m, m1)
saveRDS(meta, file='/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Li_Science/proc/meta.rds')

