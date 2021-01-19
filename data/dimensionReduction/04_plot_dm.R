#######
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/myfunc/01_function.R')
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/')
i = '2018_Li_Science'
rdir <- pdir <- paste0(i, '/dm/')
dm <- readRDS(paste0(rdir, 'diffusionMap.rds'))
rownames(dm) <- sub('.*;', '', rownames(dm))
meta <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Li_Science/proc/meta/meta.rds')
library(ggplot2)
p <- ggplot(data.frame(dm1=dm[,2],dm2=dm[,3],ct=meta$celltype),aes(x=dm1,y=dm2,col=ct)) + geom_point() + theme_classic()
ggsave(paste0(pdir, 'dm.png'), p, dpi = 200, width = 10, height = 7)

###
i = '2018_Lake_NatBiotech'
rdir <- pdir <- paste0(i, '/dm/')
dm <- readRDS(paste0(rdir, 'diffusionMap.rds'))
rownames(dm) <- sub('.*;', '', rownames(dm))
meta1 <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Lake_NatBiotech/proc/meta/CerebellarHem.rds')
meta2 <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Lake_NatBiotech/proc/meta/FrontalCortex.rds')
meta3 <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Lake_NatBiotech/proc/meta/VisualCortex.rds')
meta <- rbind(meta1, meta2, meta3)
meta[,1] <- as.character(meta[,1])
rownames(meta) <- meta[,1]
meta <- meta[rownames(dm), ]
p <- ggplot(data.frame(dm1=dm[,2],dm2=dm[,3],ct=meta$celltype),aes(x=dm1,y=dm2,col=ct)) + geom_point() + theme_classic()
ggsave(paste0(pdir, 'dm.png'), p, dpi = 200, width = 10, height = 7)

p <- ggplot(data.frame(dm1=dm[,2],dm2=dm[,3],ct=meta$celltype),aes(x=dm1,y=dm2,col=ct)) + geom_point() + theme_classic() + facet_wrap(~ct) + theme(legend.position = 'none')
ggsave(paste0(pdir, 'dm_facet.png'), p, dpi = 200, width = 30, height = 25)


###
i = '2018_Fan_CellResearch'
rdir <- pdir <- paste0(i, '/dm/')
dm <- readRDS(paste0(rdir, 'diffusionMap.rds'))
meta <- readRDS(paste0(i, '/proc/meta.rds'))
rownames(dm) <- sub('.*;', '', rownames(dm))
rownames(meta) <- meta[,1]
meta <- meta[rownames(dm), ]
p <- ggplot(data.frame(dm1=dm[,1],dm2=dm[,2],location=meta$location),aes(x=dm1,y=dm2,col=location)) + geom_point() + theme_classic()
ggsave(paste0(pdir, 'dm_location.png'), p, dpi = 200, width = 15, height = 6)
p <- ggplot(data.frame(dm2=dm[,2],dm3=dm[,3],location=meta$location),aes(x=dm2,y=dm3,col=location)) + geom_point() + theme_classic()
ggsave(paste0(pdir, 'dm_location23.png'), p, dpi = 200, width = 15, height = 6)

p <- ggplot(data.frame(dm1=dm[,1],dm2=dm[,2],time=meta$time),aes(x=dm1,y=dm2,col=time)) + geom_point() + theme_classic()
ggsave(paste0(pdir, 'dm_time.png'), p, dpi = 200, width = 15, height = 6)
p <- ggplot(data.frame(dm2=dm[,2],dm3=dm[,3],time=meta$time),aes(x=dm2,y=dm3,col=time)) + geom_point() + theme_classic()
ggsave(paste0(pdir, 'dm_time23.png'), p, dpi = 200, width = 15, height = 6)

p <- ggplot(data.frame(dm1=dm[,1],dm2=dm[,2],gender=meta$gender),aes(x=dm1,y=dm2,col=gender)) + geom_point() + theme_classic()
ggsave(paste0(pdir, 'dm_gender.png'), p, dpi = 200, width = 15, height = 6)
p <- ggplot(data.frame(dm2=dm[,2],dm3=dm[,3],gender=meta$gender),aes(x=dm2,y=dm3,col=gender)) + geom_point() + theme_classic()
ggsave(paste0(pdir, 'dm_gender23.png'), p, dpi = 200, width = 15, height = 6)

p <- ggplot(data.frame(dm1=dm[,1], dm2=dm[,2], gender=meta$gender), aes(x=dm1,y=dm2,col=gender)) + geom_point() + theme_classic()
ggsave(paste0(pdir, 'dm_gender.png'), p, dpi = 200, width = 15, height = 6)
p <- ggplot(data.frame(dm2=dm[,2], dm3=dm[,3], gender=meta$gender), aes(x=dm2,y=dm3,col=gender)) + geom_point() + theme_classic()
ggsave(paste0(pdir, 'dm_gender23.png'), p, dpi = 200, width = 15, height = 6)

p <- ggplot(data.frame(dm1=dm[,1], dm2=dm[,2], species=meta$species), aes(x=dm1,y=dm2,col=species)) + geom_point() + theme_classic()
ggsave(paste0(pdir, 'dm_species.png'), p, dpi = 200, width = 15, height = 6)
p <- ggplot(data.frame(dm2=dm[,2], dm3=dm[,3], species=meta$species), aes(x=dm2,y=dm3,col=species)) + geom_point() + theme_classic()
ggsave(paste0(pdir, 'dm_species23.png'), p, dpi = 200, width = 15, height = 6)

### not working
i = 'allen_human_fluidigm'
rdir <- pdir <- paste0(i, '/dm/')
dm <- readRDS(paste0(rdir, 'diffusionMap.rds'))
meta <- readRDS(paste0(i, '/proc/meta.rds'))
rownames(dm) <- sub('.*;', '', rownames(dm))
meta[,1] <- sub('allen', 'allen_human_fluidigm', meta[,1])
rownames(meta) <- meta[,1]
meta <- meta[rownames(dm), ]
p <- ggplot(data.frame(dm1=dm[,1], dm2=dm[,2], celltype=meta$celltype), aes(x=dm1,y=dm2,col=celltype)) + geom_point() + theme_classic()
ggsave(paste0(pdir, 'dm_celltype.png'), p, dpi = 200, width = 15, height = 6)
p <- ggplot(data.frame(dm2=dm[,2], dm3=dm[,3], celltype=meta$celltype), aes(x=dm2,y=dm3,col=celltype)) + geom_point() + theme_classic()
ggsave(paste0(pdir, 'dm_celltype23.png'), p, dpi = 200, width = 15, height = 6)

p <- ggplot(data.frame(dm1=dm[,1], dm2=dm[,2], region=meta$region), aes(x=dm1,y=dm2,col=region)) + geom_point() + theme_classic()
ggsave(paste0(pdir, 'dm_region.png'), p, dpi = 200, width = 15, height = 6)
p <- ggplot(data.frame(dm2=dm[,2], dm3=dm[,3], region=meta$region), aes(x=dm2,y=dm3,col=region)) + geom_point() + theme_classic()
ggsave(paste0(pdir, 'dm_region23.png'), p, dpi = 200, width = 15, height = 6)

p <- ggplot(data.frame(dm1=dm[,1], dm2=dm[,2], donor=meta$donor), aes(x=dm1,y=dm2,col=donor)) + geom_point() + theme_classic()
ggsave(paste0(pdir, 'dm_donor.png'), p, dpi = 200, width = 15, height = 6)
p <- ggplot(data.frame(dm2=dm[,2], dm3=dm[,3], donor=meta$donor), aes(x=dm2,y=dm3,col=donor)) + geom_point() + theme_classic()
ggsave(paste0(pdir, 'dm_donor23.png'), p, dpi = 200, width = 15, height = 6)

p <- ggplot(data.frame(dm1=dm[,1], dm2=dm[,2], layer=meta$layer), aes(x=dm1,y=dm2,col=layer)) + geom_point() + theme_classic()
ggsave(paste0(pdir, 'dm_layer.png'), p, dpi = 200, width = 15, height = 6)
p <- ggplot(data.frame(dm2=dm[,2], dm3=dm[,3], layer=meta$layer), aes(x=dm2,y=dm3,col=layer)) + geom_point() + theme_classic()
ggsave(paste0(pdir, 'dm_layer23.png'), p, dpi = 200, width = 15, height = 6)

p <- ggplot(data.frame(dm1=dm[,1], dm2=dm[,2], gender=meta$gender), aes(x=dm1,y=dm2,col=gender)) + geom_point() + theme_classic()
ggsave(paste0(pdir, 'dm_gender.png'), p, dpi = 200, width = 15, height = 6)
p <- ggplot(data.frame(dm2=dm[,2], dm3=dm[,3], gender=meta$gender), aes(x=dm2,y=dm3,col=gender)) + geom_point() + theme_classic()
ggsave(paste0(pdir, 'dm_gender23.png'), p, dpi = 200, width = 15, height = 6)


###
i <- '2017_Nowakowski_Science'## no meta

### no dm
i <- '2018_zhong_nature'
rdir <- pdir <- paste0(i, '/dm/')
dm <- readRDS(paste0(rdir, 'diffusionMap.rds')) 
meta <- readRDS(paste0(i, '/proc/meta.rds'))
rownames(dm) <- sub('.*;', '', rownames(dm))
rownames(meta) <- meta[,1]
meta <- meta[rownames(dm), ]
meta <- as.data.frame(meta)
p <- ggplot(data.frame(dm1=dm[,1], dm2=dm[,2], Gender=meta$Gender), aes(x=dm1,y=dm2,col=Gender)) + geom_point() + theme_classic()
ggsave(paste0(pdir, 'dm_Gender.png'), p, dpi = 200, width = 12, height = 6)
p <- ggplot(data.frame(dm2=dm[,2], dm3=dm[,3], Gender=meta$Gender), aes(x=dm2,y=dm3,col=Gender)) + geom_point() + theme_classic()
ggsave(paste0(pdir, 'dm_Gender23.png'), p, dpi = 200, width = 12, height = 6)

p <- ggplot(data.frame(dm1=dm[,1], dm2=dm[,2], Location=meta$Location), aes(x=dm1,y=dm2,col=Location)) + geom_point() + theme_classic()
ggsave(paste0(pdir, 'dm_Location.png'), p, dpi = 200, width = 12, height = 6)
p <- ggplot(data.frame(dm2=dm[,2], dm3=dm[,3], Location=meta$Location), aes(x=dm2,y=dm3,col=Location)) + geom_point() + theme_classic()
ggsave(paste0(pdir, 'dm_Location23.png'), p, dpi = 200, width = 12, height = 6)

p <- ggplot(data.frame(dm1=dm[,1], dm2=dm[,2], Time=meta$Time), aes(x=dm1,y=dm2,col=Time)) + geom_point() + theme_classic()
ggsave(paste0(pdir, 'dm_Time.png'), p, dpi = 200, width = 12, height = 6)
p <- ggplot(data.frame(dm2=dm[,2], dm3=dm[,3], Time=meta$Time), aes(x=dm2,y=dm3,col=Time)) + geom_point() + theme_classic()
ggsave(paste0(pdir, 'dm_Time23.png'), p, dpi = 200, width = 12, height = 6)

###
i <- '2015_Darmanis_PNAS'
rdir <- pdir <- paste0(i, '/dm/')
dm <- readRDS(paste0(rdir, 'diffusionMap.rds')) ## no meta
meta <- readRDS(paste0(i, '/proc/meta.rds'))
rownames(dm) <- sub('.*;', '', rownames(dm))
rownames(meta) <- meta[,1]
meta <- meta[rownames(dm), ]
meta <- as.data.frame(meta)
p <- ggplot(data.frame(dm1=dm[,1], dm2=dm[,2], celltype=meta$celltype), aes(x=dm1,y=dm2,col=celltype)) + geom_point() + theme_classic()
ggsave(paste0(pdir, 'dm_celltype.png'), p, dpi = 200, width = 12, height = 6)
p <- ggplot(data.frame(dm2=dm[,2], dm3=dm[,3], celltype=meta$celltype), aes(x=dm2,y=dm3,col=celltype)) + geom_point() + theme_classic()
ggsave(paste0(pdir, 'dm_celltype23.png'), p, dpi = 200, width = 12, height = 6)

p <- ggplot(data.frame(dm1=dm[,1], dm2=dm[,2], time=meta$time), aes(x=dm1,y=dm2,col=time)) + geom_point() + theme_classic()
ggsave(paste0(pdir, 'dm_time.png'), p, dpi = 200, width = 12, height = 6)
p <- ggplot(data.frame(dm2=dm[,2], dm3=dm[,3], time=meta$time), aes(x=dm2,y=dm3,col=time)) + geom_point() + theme_classic()
ggsave(paste0(pdir, 'dm_time23.png'), p, dpi = 200, width = 12, height = 6)

p <- ggplot(data.frame(dm1=dm[,1], dm2=dm[,2], species=meta$species), aes(x=dm1,y=dm2,col=species)) + geom_point() + theme_classic()
ggsave(paste0(pdir, 'dm_species.png'), p, dpi = 200, width = 12, height = 6)
p <- ggplot(data.frame(dm2=dm[,2], dm3=dm[,3], species=meta$species), aes(x=dm2,y=dm3,col=species)) + geom_point() + theme_classic()
ggsave(paste0(pdir, 'dm_species23.png'), p, dpi = 200, width = 12, height = 6)

p <- ggplot(data.frame(dm1=dm[,1], dm2=dm[,2], location=meta$location), aes(x=dm1,y=dm2,col=location)) + geom_point() + theme_classic()
ggsave(paste0(pdir, 'dm_location.png'), p, dpi = 200, width = 12, height = 6)
p <- ggplot(data.frame(dm2=dm[,2], dm3=dm[,3], location=meta$location), aes(x=dm2,y=dm3,col=location)) + geom_point() + theme_classic()
ggsave(paste0(pdir, 'dm_location23.png'), p, dpi = 200, width = 12, height = 6)


###
i <- '2014_Pollen_NatBiotech'
rdir <- pdir <- paste0(i, '/dm/')
dm <- readRDS(paste0(rdir, 'diffusionMap.rds')) ## no meta
meta <- readRDS(paste0(i, '/proc/meta.rds'))
rownames(dm) <- sub('.*;', '', rownames(dm))
rownames(meta) <- meta[,1]
meta <- meta[rownames(dm), ]
meta <- as.data.frame(meta)
p <- ggplot(data.frame(dm1=dm[,1], dm2=dm[,2], celltype=meta$celltype), aes(x=dm1,y=dm2,col=celltype)) + geom_point() + theme_classic()
ggsave(paste0(pdir, 'dm_celltype.png'), p, dpi = 200, width = 12, height = 6)
p <- ggplot(data.frame(dm2=dm[,2], dm3=dm[,3], celltype=meta$celltype), aes(x=dm2,y=dm3,col=celltype)) + geom_point() + theme_classic()
ggsave(paste0(pdir, 'dm_celltype23.png'), p, dpi = 200, width = 12, height = 6)

###
i <- '2020_Guo_Nature'
rdir <- pdir <- paste0(i, '/dm/')
dm <- readRDS(paste0(rdir, 'diffusionMap.rds')) ## no meta
meta <- readRDS(paste0(i, '/meta/meta_sub.rds'))
rownames(dm) <- sub('.*;', '', rownames(dm))
rownames(meta) <- meta[,1]
meta <- meta[rownames(dm), ]
meta <- as.data.frame(meta)
p <- ggplot(data.frame(dm1=dm[,1], dm2=dm[,2], celltype=meta$celltype), aes(x=dm1,y=dm2,col=celltype)) + geom_point() + theme_classic()
ggsave(paste0(pdir, 'dm_celltype.png'), p, dpi = 200, width = 12, height = 6)
p <- ggplot(data.frame(dm2=dm[,2], dm3=dm[,3], celltype=meta$celltype), aes(x=dm2,y=dm3,col=celltype)) + geom_point() + theme_classic()
ggsave(paste0(pdir, 'dm_celltype23.png'), p, dpi = 200, width = 12, height = 6)
p <- ggplot(data.frame(dm1=dm[,2], dm3=dm[,3], celltype=meta$celltype), aes(x=dm1,y=dm3,col=celltype)) + geom_point() + theme_classic()
ggsave(paste0(pdir, 'dm_celltype13.png'), p, dpi = 200, width = 12, height = 6)


p <- ggplot(data.frame(dm1=dm[,1], dm2=dm[,2], celltype=meta$celltype), aes(x=dm1,y=dm2,col=celltype)) + geom_point() + theme_classic() + facet_wrap(~celltype)
ggsave(paste0(pdir, 'dm_celltype_facet.png'), p, dpi = 200, width = 20, height = 15)
p <- ggplot(data.frame(dm2=dm[,2], dm3=dm[,3], celltype=meta$celltype), aes(x=dm2,y=dm3,col=celltype)) + geom_point() + theme_classic()+ facet_wrap(~celltype)
ggsave(paste0(pdir, 'dm_celltype23_facet.png'), p, dpi = 200, width = 20, height = 15)
p <- ggplot(data.frame(dm1=dm[,2], dm3=dm[,3], celltype=meta$celltype), aes(x=dm1,y=dm3,col=celltype)) + geom_point() + theme_classic()+ facet_wrap(~celltype)
ggsave(paste0(pdir, 'dm_celltype13_facet.png'), p, dpi = 200, width = 20, height = 15)

p <- ggplot(data.frame(dm1=dm[,1], dm2=dm[,2], time=meta$time), aes(x=dm1,y=dm2,col=time)) + geom_point() + theme_classic()
ggsave(paste0(pdir, 'dm_time.png'), p, dpi = 200, width = 12, height = 6)
p <- ggplot(data.frame(dm2=dm[,2], dm3=dm[,3], time=meta$time), aes(x=dm2,y=dm3,col=time)) + geom_point() + theme_classic()
ggsave(paste0(pdir, 'dm_time23.png'), p, dpi = 200, width = 12, height = 6)

p <- ggplot(data.frame(dm1=dm[,1], dm2=dm[,2], donor=meta$donor), aes(x=dm1,y=dm2,col=donor)) + geom_point() + theme_classic()
ggsave(paste0(pdir, 'dm_donor.png'), p, dpi = 200, width = 12, height = 6)
p <- ggplot(data.frame(dm2=dm[,2], dm3=dm[,3], donor=meta$donor), aes(x=dm2,y=dm3,col=donor)) + geom_point() + theme_classic()
ggsave(paste0(pdir, 'dm_donor23.png'), p, dpi = 200, width = 12, height = 6)

p <- ggplot(data.frame(dm1=dm[,1], dm2=dm[,2], study=meta$study), aes(x=dm1,y=dm2,col=study)) + geom_point() + theme_classic()
ggsave(paste0(pdir, 'dm_study.png'), p, dpi = 200, width = 12, height = 6)
p <- ggplot(data.frame(dm2=dm[,2], dm3=dm[,3], study=meta$study), aes(x=dm2,y=dm3,col=study)) + geom_point() + theme_classic()
ggsave(paste0(pdir, 'dm_study23.png'), p, dpi = 200, width = 12, height = 6)
 
library("gg3D")
p <- ggplot(data.frame(dm1=dm[,1], dm2=dm[,2], dm3=dm[,3], celltype=meta$celltype), aes(x=dm1, y=dm2, z=dm3, color=celltype)) + 
  theme_void() +
  axes_3D() +
  stat_3D()
ggsave(paste0(pdir, 'dm_celltype_3D.png'), p, dpi = 200, width = 12, height = 6)

## allen_human_10x         




