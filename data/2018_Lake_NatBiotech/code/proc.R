library(GEOsearch)
library(data.table)
m <- SampleDetail('GSE97930')
brodmann <- sub('brodmann area: ','',sapply(m[,'Characteristic'],function(i) grep('brodmann area:',strsplit(i,'; ')[[1]],value = T),USE.NAMES = F))
age <- sub('age: ','',sapply(m[,'Characteristic'],function(i) grep('age:',strsplit(i,'; ')[[1]],value = T),USE.NAMES = F))
sex <- sub('Sex: ','',sapply(m[,'Characteristic'],function(i) grep('Sex:',strsplit(i,'; ')[[1]],value = T),USE.NAMES = F))
race <- sub('race: ','',sapply(m[,'Characteristic'],function(i) grep('race:',strsplit(i,'; ')[[1]],value = T),USE.NAMES = F))
region <- sub('region: ','',sapply(m[,'Characteristic'],function(i) grep('region:',strsplit(i,'; ')[[1]],value = T),USE.NAMES = F))
brodmann[brodmann=='Cerebellar Hemisphere'] <- NA
region[region=='Lateral'] <- 'Lateral Cerebellar Hemisphere'
region <- sub('brain ','',region)

sampm <- data.frame(title=tolower(m$Title),brodmann,age,sex,race,region,stringsAsFactors = F)
sampm <- rbind(sampm,data.frame(title='Cer',brodmann=NA,age=NA,sex=NA,race=NA,region='Lateral Cerebellar Hemisphere',stringsAsFactors = F))


conv <- c('Ex'='excitatory neurons','In'='inhibitory neurons','Gran'= 'granule cells','Purk1'= 'Purkinje neurons','Purk2'= 'Purkinje neurons', 'End'='endothelial cells', 'Ast'= 'astrocytes', 'Oli'='oligodendrocytes', 'OPC'= 'oligodendrocytes precursor cells', 'Mic'='microglia','Per'='pericytes')

for (type in c('CerebellarHem','FrontalCortex','VisualCortex')) {
  d=fread(paste0('/scratch/users/whou10@jhu.edu/Wenpin/brain/data/2018_Lake_NatBiotech/raw/GSE97930_',type,'_snDrop-seq_UMI_Count_Matrix_08-01-2017.txt'),data.table=F)
  rownames(d) = d[,1]
  d = d[,-1]
  d = as.matrix(d)
  ct <- sub('_.*','',colnames(d))
  cct <- ct
  cct[cct %in% names(conv)] <- conv[cct[cct %in% names(conv)]]
  cct[grep('Ex[1-9].*',ct)] <- 'excitatory neurons'
  cct[grep('In[1-9].*',ct)] <- 'inhibitory neurons'
  
  c2 <- sapply(colnames(d),function(i) strsplit(i,'_')[[1]][2],USE.NAMES = F)
  cn <- paste0('2018_Lake_NatBiotech:',type,':cell_',1:ncol(d))
  colnames(d) <- cn
  meta <- data.frame(cell=cn,celltype=cct,sampm[match(c2,sampm[,1]),],species='human')
  saveRDS(d,paste0("/scratch/users/whou10@jhu.edu/Wenpin/brain/data/2018_Lake_NatBiotech/proc/count/",type,".rds"))
  saveRDS(meta,paste0("/scratch/users/whou10@jhu.edu/Wenpin/brain/data/2018_Lake_NatBiotech/proc/meta/",type,".rds"))
}


