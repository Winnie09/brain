setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/brain')
pos = readRDS('./pseudotime/plot/plot/pseudotime_order.rds')
pairwiseKLD <- function(group1, group2){
  min = range(d[,'pseudotime'])[1]
  max = range(d[,'pseudotime'])[2]
  as.vector(sapply(group1, function(p1){
        sapply(group2, function(p2){
          a = density(d[d$Sample.ID==p1,'pseudotime'], from = min, to = max)
          b = density(d[d$Sample.ID==p2,'pseudotime'], from = min, to = max)
          KLD(a$y,b$y)$mean.sum.KLD
        })
  }))
}
groupKLD <- function(group1, group2){
  v1 = pairwiseKLD(group1, group2)
  v2 = as.vector(sapply(1:1000, function(i){
    if (!i%%10) print(i)
    set.seed(i)
    ap = c(group1, group2)
    gp1 = sample(ap, length(group1))
    gp2 = sample(ap, length(group2))
    mean(pairwiseKLD(gp1, gp2))  ## use mean only in permuation test
  }))
  ## KLD
  # KLD(density(v1, from = min(v1,v2), to = max(v1,v2))$y, density(v2, from = min(v1,v2), to = max(v1,v2))$y)$mean.sum.KLD
  ## t test
  # t.test(v1,v2)$p.value
  ## permutation test
  sum(v2 > mean(v1))/length(v2)
}
featureKLD <- function(feature){
  aval = unique(as.character(d[,feature]))
  m = sapply(aval, function(val1){
    sapply(aval, function(val2){
        group1 = unique(d[as.character(d[,feature])==val1,'Sample.ID'])
        group2 = unique(d[as.character(d[,feature])==val2,'Sample.ID'])
        groupKLD(group1,group2)      
    })
  })
  m
}

alltraj = names(pos)

library(gplots)
for (traj in setdiff(alltraj,'6_4')){
  print(traj)
  plotdir = paste0('./pseudotime_variability/plot/', traj,'/')
  dir.create(plotdir, recursive = T, showWarnings = F)
  d = pos[[traj]]
  allfeat = setdiff(colnames(d)[5:(ncol(d)-1)],c('Pathology','Tumor.Grade','Age'))
  pseudotime <- 'pseudotime'
  for (feature in allfeat){
    print(feature)
    # feature = 'Location'
    pdf(paste0(plotdir,feature,'_time_distribution.pdf'))
    print(ggplot(d,aes_string(x=pseudotime,fill=feature),alpha=0.5) + geom_histogram() + theme_classic() + facet_wrap(~patient,scales = 'free_y',ncol=1,strip.position="right") + theme(strip.text.y = element_text(angle = 0)) + 
            theme(axis.text=element_text(size=5))+
            scale_fill_brewer(palette='Set2'))
    dev.off()
    if (length(unique(as.character(d[,feature])))>1){
      x=featureKLD(feature)
      pdf(paste0(plotdir,feature,'.pdf'),width=6,height=6)
      print(heatmap.2(x))
      dev.off()  
    }
    
  }
}


