af <- list.files(paste0('/home-4/zji4@jhu.edu/scratch/cnvsel/cnv/'),pattern = '.rds')
tmp <- sapply(af,function(f) {
  print(f)
  d <- readRDS(paste0('/home-4/zji4@jhu.edu/scratch/cnvsel/cnv/',f))
  print(sapply(d,function(i) sapply(i,length)))
})
