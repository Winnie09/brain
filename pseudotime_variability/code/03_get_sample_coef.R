order = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/pseudotime_variability/result/6_10_3/order.rds')
order = order[grepl('GBM', order$sample_name),]
ap = gsub('_.*','', order$sample_name)

gbm <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/GBM/data/proc/tumor/matrix/saver.rds')
gbm = gbm[, order$sample_name]

source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R')
res <- sapply(unique(ap), function(s){
  trainData = gbm[, grepl(s, colnames(gbm))]
  trainX = order[grepl(s, order$sample_name),'Pseudotime']
  para = get_spline_coefficient(trainData=trainData, trainX=trainX, fit.min=min(order$Pseudotime), fit.max=max(order$Pseudotime))
  saveRDS(para, paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/pseudotime_variability/result/6_10_3/spline_coef_', s, '.rds'))
  return(0)
})

