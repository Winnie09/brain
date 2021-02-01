#######
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/myfunc/01_function.R')
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/')
i = '2018_Li_Science'
rdir <- pdir <- paste0(i, '/umap/')
umap <- readRDS(paste0(rdir, 'umap.rds'))
rownames(umap) <- sub('.*;', '', rownames(umap))
meta <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Li_Science/proc/meta/meta.rds')
library(ggplot2)
p <- ggplot(data.frame(umap1=umap[,1],umap2=umap[,2],ct=meta$celltype),aes(x=umap1,y=umap2,col=ct)) + geom_point() + theme_classic()
ggsave(paste0(pdir, 'umap.png'), p, dpi = 200, width = 10, height = 7)
p <- ggplot(data.frame(umap1=umap[,1],umap2=umap[,2],ct=meta$celltype),aes(x=umap1,y=umap2,col=ct)) + geom_point() + theme_classic() + facet_wrap(~ct) + theme(legend.position = 'none')
ggsave(paste0(pdir, 'umap_facet.png'), p, dpi = 200, width = 30, height = 25)

###
i = '2018_Lake_NatBiotech'
rdir <- pdir <- paste0(i, '/umap/')
umap <- readRDS(paste0(rdir, 'umap.rds'))
rownames(umap) <- sub('.*;', '', rownames(umap))
meta1 <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Lake_NatBiotech/proc/meta/CerebellarHem.rds')
meta2 <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Lake_NatBiotech/proc/meta/FrontalCortex.rds')
meta3 <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/data/2018_Lake_NatBiotech/proc/meta/VisualCortex.rds')
meta <- rbind(meta1, meta2, meta3)
meta[,1] <- as.character(meta[,1])
rownames(meta) <- meta[,1]
meta <- meta[rownames(umap), ]
p <- ggplot(data.frame(umap1=umap[,1],umap2=umap[,2],ct=meta$celltype),aes(x=umap1,y=umap2,col=ct)) + geom_point() + theme_classic()
ggsave(paste0(pdir, 'umap.png'), p, dpi = 200, width = 10, height = 7)
p <- ggplot(data.frame(umap1=umap[,1],umap2=umap[,2],ct=meta$celltype),aes(x=umap1,y=umap2,col=ct)) + geom_point() + theme_classic() + facet_wrap(~ct) + theme(legend.position = 'none')
ggsave(paste0(pdir, 'umap_facet.png'), p, dpi = 200, width = 30, height = 25)


###
i = '2018_Fan_CellResearch'
rdir <- pdir <- paste0(i, '/umap/')
umap <- readRDS(paste0(rdir, 'umap.rds'))
meta <- readRDS(paste0(i, '/proc/meta.rds'))
rownames(umap) <- sub('.*;', '', rownames(umap))
rownames(meta) <- meta[,1]
meta <- meta[rownames(umap), ]
p <- ggplot(data.frame(umap1=umap[,1],umap2=umap[,2],location=meta$location),aes(x=umap1,y=umap2,col=location)) + geom_point() + theme_classic()
ggsave(paste0(pdir, 'umap_location.png'), p, dpi = 200, width = 15, height = 6)
p <- ggplot(data.frame(umap1=umap[,1],umap2=umap[,2],location=meta$location),aes(x=umap1,y=umap2,col=location)) + geom_point() + theme_classic() + facet_wrap(~location) + theme(legend.position = 'none')
ggsave(paste0(pdir, 'umap_location_facet.png'), p, dpi = 200, width = 30, height = 25)


p <- ggplot(data.frame(umap1=umap[,1],umap2=umap[,2],time=meta$time),aes(x=umap1,y=umap2,col=time)) + geom_point() + theme_classic()
ggsave(paste0(pdir, 'umap_time.png'), p, dpi = 200, width = 15, height = 6)
p <- ggplot(data.frame(umap1=umap[,1],umap2=umap[,2],time=meta$time),aes(x=umap1,y=umap2,col=time)) + geom_point() + theme_classic() +facet_wrap(~time) + theme(legend.position = 'none')
ggsave(paste0(pdir, 'umap_time_facet.png'), p, dpi = 200, width = 30, height = 25)


p <- ggplot(data.frame(umap1=umap[,1],umap2=umap[,2],gender=meta$gender),aes(x=umap1,y=umap2,col=gender)) + geom_point() + theme_classic()
ggsave(paste0(pdir, 'umap_gender.png'), p, dpi = 200, width = 15, height = 6)
p <- ggplot(data.frame(umap1=umap[,1],umap2=umap[,2],gender=meta$gender),aes(x=umap1,y=umap2,col=gender)) + geom_point() + theme_classic() + facet_wrap(~gender) + theme(legend.position = 'none')
ggsave(paste0(pdir, 'umap_gender_facet.png'), p, dpi = 200, width = 30, height = 25)



p <- ggplot(data.frame(umap1=umap[,1], umap2=umap[,2], species=meta$species), aes(x=umap1,y=umap2,col=species)) + geom_point() + theme_classic()
ggsave(paste0(pdir, 'umap_species.png'), p, dpi = 200, width = 15, height = 6)
p <- ggplot(data.frame(umap1=umap[,1], umap2=umap[,2], species=meta$species), aes(x=umap1,y=umap2,col=species)) + geom_point() + theme_classic() + facet_wrap(~species) + theme(legend.position = 'none')
ggsave(paste0(pdir, 'umap_species_facet.png'), p, dpi = 200, width = 30, height = 25)

### not working
i = 'allen_human_fluidigm'
rdir <- pdir <- paste0(i, '/umap/')
umap <- readRDS(paste0(rdir, 'umap.rds'))
meta <- readRDS(paste0(i, '/proc/meta.rds'))
rownames(umap) <- sub('.*;', '', rownames(umap))
meta[,1] <- sub('allen', 'allen_human_fluidigm', meta[,1])
rownames(meta) <- meta[,1]
meta <- meta[rownames(umap), ]
p <- ggplot(data.frame(umap1=umap[,1], umap2=umap[,2], celltype=meta$celltype), aes(x=umap1,y=umap2,col=celltype)) + geom_point() + theme_classic()
ggsave(paste0(pdir, 'umap_celltype.png'), p, dpi = 200, width = 15, height = 6)
p <- ggplot(data.frame(umap1=umap[,1], umap2=umap[,2], celltype=meta$celltype), aes(x=umap1,y=umap2,col=celltype)) + geom_point() + theme_classic() + facet_wrap(~celltype) + theme(legend.position = 'none')
ggsave(paste0(pdir, 'umap_celltype_facet.png'), p, dpi = 200, width = 30, height = 25)


p <- ggplot(data.frame(umap1=umap[,1], umap2=umap[,2], region=meta$region), aes(x=umap1,y=umap2,col=region)) + geom_point() + theme_classic()
ggsave(paste0(pdir, 'umap_region.png'), p, dpi = 200, width = 15, height = 6)
p <- ggplot(data.frame(umap1=umap[,1], umap2=umap[,2], region=meta$region), aes(x=umap1,y=umap2,col=region)) + geom_point() + theme_classic() + facet_wrap(~region) + theme(legend.position = 'none')
ggsave(paste0(pdir, 'umap_region_facet.png'), p, dpi = 200, width = 30, height = 25)

p <- ggplot(data.frame(umap1=umap[,1], umap2=umap[,2], donor=meta$donor), aes(x=umap1,y=umap2,col=donor)) + geom_point() + theme_classic()
ggsave(paste0(pdir, 'umap_donor.png'), p, dpi = 200, width = 15, height = 6)
p <- ggplot(data.frame(umap1=umap[,1], umap2=umap[,2], donor=meta$donor), aes(x=umap1,y=umap2,col=donor)) + geom_point() + theme_classic() + facet_wrap(~donor) + theme(legend.position = 'none')
ggsave(paste0(pdir, 'umap_donor_facet.png'), p, dpi = 200, width = 30, height = 25)


p <- ggplot(data.frame(umap1=umap[,1], umap2=umap[,2], layer=meta$layer), aes(x=umap1,y=umap2,col=layer)) + geom_point() + theme_classic()
ggsave(paste0(pdir, 'umap_layer.png'), p, dpi = 200, width = 15, height = 6)
p <- ggplot(data.frame(umap1=umap[,1], umap2=umap[,2], layer=meta$layer), aes(x=umap1,y=umap2,col=layer)) + geom_point() + theme_classic() + facet_wrap(~layer) + theme(legend.position = 'none')
ggsave(paste0(pdir, 'umap_layer_facet.png'), p, dpi = 200, width = 30, height = 25)

p <- ggplot(data.frame(umap1=umap[,1], umap2=umap[,2], gender=meta$gender), aes(x=umap1,y=umap2,col=gender)) + geom_point() + theme_classic()
ggsave(paste0(pdir, 'umap_gender.png'), p, dpi = 200, width = 15, height = 6)
p <- ggplot(data.frame(umap1=umap[,1], umap2=umap[,2], gender=meta$gender), aes(x=umap1,y=umap2,col=gender)) + geom_point() + theme_classic() + theme_classic() + facet_wrap(~gender) + theme(legend.position = 'none')
ggsave(paste0(pdir, 'umap_gender_facet.png'), p, dpi = 200, width = 30, height = 25)


###
i <- '2017_Nowakowski_Science'## no meta

### no umap
i <- '2018_zhong_nature'
rdir <- pdir <- paste0(i, '/umap/')
umap <- readRDS(paste0(rdir, 'umap.rds')) 
meta <- readRDS(paste0(i, '/proc/meta.rds'))
rownames(umap) <- sub('.*;', '', rownames(umap))
rownames(meta) <- meta[,1]
meta <- meta[rownames(umap), ]
meta <- as.data.frame(meta)
p <- ggplot(data.frame(umap1=umap[,1], umap2=umap[,2], Gender=meta$Gender), aes(x=umap1,y=umap2,col=Gender)) + geom_point() + theme_classic()
ggsave(paste0(pdir, 'umap_Gender.png'), p, dpi = 200, width = 12, height = 6)
p <- ggplot(data.frame(umap1=umap[,1], umap2=umap[,2], Gender=meta$Gender), aes(x=umap1,y=umap2,col=Gender)) + geom_point() + theme_classic() + facet_wrap(~Gender) + theme(legend.position = 'none')
ggsave(paste0(pdir, 'umap_Gender_facet.png'), p, dpi = 200, width = 30, height = 25)

p <- ggplot(data.frame(umap1=umap[,1], umap2=umap[,2], Location=meta$Location), aes(x=umap1,y=umap2,col=Location)) + geom_point() + theme_classic()
ggsave(paste0(pdir, 'umap_Location.png'), p, dpi = 200, width = 12, height = 6)
p <- ggplot(data.frame(umap1=umap[,1], umap2=umap[,2], Location=meta$Location), aes(x=umap1,y=umap2,col=Location)) + geom_point() + theme_classic() + facet_wrap(~Location) + theme(legend.position = 'none')
ggsave(paste0(pdir, 'umap_Location_facet.png'), p, dpi = 200, width = 30, height = 25)

p <- ggplot(data.frame(umap1=umap[,1], umap2=umap[,2], Time=meta$Time), aes(x=umap1,y=umap2,col=Time)) + geom_point() + theme_classic()
ggsave(paste0(pdir, 'umap_Time.png'), p, dpi = 200, width = 12, height = 6)
p <- ggplot(data.frame(umap1=umap[,1], umap2=umap[,2], Time=meta$Time), aes(x=umap1,y=umap2,col=Time)) + geom_point() + theme_classic() + facet_wrap(~Time) + theme(legend.position = 'none')
ggsave(paste0(pdir, 'umap_Time_facet.png'), p, dpi = 200, width = 30, height = 25)

###
i <- '2015_Darmanis_PNAS'
rdir <- pdir <- paste0(i, '/umap/')
umap <- readRDS(paste0(rdir, 'umap.rds')) ## no meta
meta <- readRDS(paste0(i, '/proc/meta.rds'))
rownames(umap) <- sub('.*;', '', rownames(umap))
rownames(meta) <- meta[,1]
meta <- meta[rownames(umap), ]
meta <- as.data.frame(meta)
p <- ggplot(data.frame(umap1=umap[,1], umap2=umap[,2], celltype=meta$celltype), aes(x=umap1,y=umap2,col=celltype)) + geom_point() + theme_classic()
ggsave(paste0(pdir, 'umap_celltype.png'), p, dpi = 200, width = 12, height = 6)
p <- ggplot(data.frame(umap1=umap[,1], umap2=umap[,2], celltype=meta$celltype), aes(x=umap1,y=umap2,col=celltype)) + geom_point() + theme_classic() + facet_wrap(~celltype) + theme(legend.position = 'none')
ggsave(paste0(pdir, 'umap_celltype_facet.png'), p, dpi = 200, width = 30, height = 25)


p <- ggplot(data.frame(umap1=umap[,1], umap2=umap[,2], time=meta$time), aes(x=umap1,y=umap2,col=time)) + geom_point() + theme_classic()
ggsave(paste0(pdir, 'umap_time.png'), p, dpi = 200, width = 12, height = 6)
p <- ggplot(data.frame(umap1=umap[,1], umap2=umap[,2], time=meta$time), aes(x=umap1,y=umap2,col=time)) + geom_point() + theme_classic() + facet_wrap(~time) + theme(legend.position = 'none')
ggsave(paste0(pdir, 'umap_time_facet.png'), p, dpi = 200, width = 30, height = 25)


p <- ggplot(data.frame(umap1=umap[,1], umap2=umap[,2], species=meta$species), aes(x=umap1,y=umap2,col=species)) + geom_point() + theme_classic()
ggsave(paste0(pdir, 'umap_species.png'), p, dpi = 200, width = 12, height = 6)
p <- ggplot(data.frame(umap1=umap[,1], umap2=umap[,2], species=meta$species), aes(x=umap1,y=umap2,col=species)) + geom_point() + theme_classic() + facet_wrap(~species) + theme(legend.position = 'none')
ggsave(paste0(pdir, 'umap_species_facet.png'), p, dpi = 200, width = 30, height = 25)

p <- ggplot(data.frame(umap1=umap[,1], umap2=umap[,2], location=meta$location), aes(x=umap1,y=umap2,col=location)) + geom_point() + theme_classic()
ggsave(paste0(pdir, 'umap_location.png'), p, dpi = 200, width = 12, height = 6)
p <- ggplot(data.frame(umap1=umap[,1], umap2=umap[,2], location=meta$location), aes(x=umap1,y=umap2,col=location)) + geom_point() + theme_classic() + facet_wrap(~location) + theme(legend.position = 'none')
ggsave(paste0(pdir, 'umap_location_facet.png'), p, dpi = 200, width = 30, height = 25)



###
i <- '2014_Pollen_NatBiotech'
rdir <- pdir <- paste0(i, '/umap/')
umap <- readRDS(paste0(rdir, 'umap.rds')) ## no meta
meta <- readRDS(paste0(i, '/proc/meta.rds'))
rownames(umap) <- sub('.*;', '', rownames(umap))
rownames(meta) <- meta[,1]
meta <- meta[rownames(umap), ]
meta <- as.data.frame(meta)
p <- ggplot(data.frame(umap1=umap[,1], umap2=umap[,2], celltype=meta$celltype), aes(x=umap1,y=umap2,col=celltype)) + geom_point() + theme_classic()
ggsave(paste0(pdir, 'umap_celltype.png'), p, dpi = 200, width = 12, height = 6)
p <- ggplot(data.frame(umap1=umap[,1], umap2=umap[,2], celltype=meta$celltype), aes(x=umap1,y=umap2,col=celltype)) + geom_point() + theme_classic() + facet_wrap(~celltype) + theme(legend.position = 'none')
ggsave(paste0(pdir, 'umap_celltype_facet.png'), p, dpi = 200, width = 30, height = 25)

###
i <- '2020_Guo_Nature'
rdir <- pdir <- paste0(i, '/umap/')
umap <- readRDS(paste0(rdir, 'umap.rds')) ## no meta
meta <- readRDS(paste0(i, '/meta/meta_sub.rds'))
rownames(umap) <- sub('.*;', '', rownames(umap))
rownames(meta) <- meta[,1]
meta <- meta[rownames(umap), ]
meta <- as.data.frame(meta)
p <- ggplot(data.frame(umap1=umap[,1], umap2=umap[,2], celltype=meta$celltype), aes(x=umap1,y=umap2,col=celltype)) + geom_point() + theme_classic()
ggsave(paste0(pdir, 'umap_celltype.png'), p, dpi = 200, width = 12, height = 6)
p <- ggplot(data.frame(umap1=umap[,1], umap2=umap[,2], celltype=meta$celltype), aes(x=umap1,y=umap2,col=celltype)) + geom_point() + theme_classic() + facet_wrap(~celltype) + theme(legend.position = 'none')
ggsave(paste0(pdir, 'umap_celltype_facet.png'), p, dpi = 200, width = 30, height = 25)

p <- ggplot(data.frame(umap1=umap[,1], umap2=umap[,2], time=meta$time), aes(x=umap1,y=umap2,col=time)) + geom_point() + theme_classic()
ggsave(paste0(pdir, 'umap_time.png'), p, dpi = 200, width = 12, height = 6)
p <- ggplot(data.frame(umap1=umap[,1], umap2=umap[,2], time=meta$time), aes(x=umap1,y=umap2,col=time)) + geom_point() + theme_classic()+ facet_wrap(~time) + theme(legend.position = 'none')
ggsave(paste0(pdir, 'umap_time_facet.png'), p, dpi = 200, width = 12, height = 6)


p <- ggplot(data.frame(umap1=umap[,1], umap2=umap[,2], donor=meta$donor), aes(x=umap1,y=umap2,col=donor)) + geom_point() + theme_classic()
ggsave(paste0(pdir, 'umap_donor.png'), p, dpi = 200, width = 12, height = 6)
p <- ggplot(data.frame(umap1=umap[,1], umap2=umap[,2], donor=meta$donor), aes(x=umap1,y=umap2,col=donor)) + geom_point() + theme_classic()+ facet_wrap(~donor) + theme(legend.position = 'none')
ggsave(paste0(pdir, 'umap_donor_facet.png'), p, dpi = 200, width = 12, height = 6)


p <- ggplot(data.frame(umap1=umap[,1], umap2=umap[,2], study=meta$study), aes(x=umap1,y=umap2,col=study)) + geom_point() + theme_classic()
ggsave(paste0(pdir, 'umap_study.png'), p, dpi = 200, width = 12, height = 6)
p <- ggplot(data.frame(umap1=umap[,1], umap2=umap[,2], study=meta$study), aes(x=umap1,y=umap2,col=study)) + geom_point() + theme_classic()+ facet_wrap(~study) + theme(legend.position = 'none')
ggsave(paste0(pdir, 'umap_study_facet.png'), p, dpi = 200, width = 12, height = 6)

## allen_human_10x         

