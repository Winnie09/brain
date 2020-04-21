par(mfrow=c(2,2))
for (seed in 1:4){
  set.seed(seed)
  trainData = matrix(rnorm(20000), ncol=100)
  trainData[sample(1:200,50*(seed-1)), sample(1:100,20*(seed-1))] <- trainData[sample(1:200,50*(seed-1)), sample(1:100,20*(seed-1))] -10
  # trainData = scale(trainData)
  time = seq(1,100)
  trainX = time
  num.base = 10
  knots = seq(min(time),max(time),length.out=num.base+2)[2:(num.base+1)]
  library(splines)
  base = cbind(1,bs(trainX,knots = knots))
  colidx = NULL
  for (ii in seq(2,ncol(base))){
    if (length(unique(base[,ii])) == 1) colidx = c(colidx, ii)
  }
  if (length(colidx)) base = base[,-colidx]
  coef = t(chol2inv(chol(crossprod(base))) %*% t(base) %*% t(trainData))
  # pred = predict(lm(trainData[1,]~base-1))
  pred = base %*% coef[1,]
  
  pd = data.frame(x=time,y=trainData[1,],pred=pred)
  plot(pd$y~pd$x,pch=20,cex=1)
  lines(pred,col='red',lwd=2)
}

#####
par(mfrow=c(2,2))
for (seed in 1:4){
  set.seed(seed)
  trainData = matrix(rnorm(20000), ncol=100)
  
  time = seq(1,100)
  trainX = time
  num.base = 10
  knots = seq(min(time),max(time),length.out=num.base+2)[2:(num.base+1)]
  library(splines)
  base = cbind(1,bs(trainX,knots = knots))
  colidx = NULL
  for (ii in seq(2,ncol(base))){
    if (length(unique(base[,ii])) == 1) colidx = c(colidx, ii)
  }
  if (length(colidx)) base = base[,-colidx]
  coef = t(chol2inv(chol(crossprod(base))) %*% t(base) %*% t(trainData))
  pred = base %*% coef[1,]
    
  pd = data.frame(x=time,y=trainData[1,],pred=pred)
  plot(pd$y~pd$x,pch=20,cex=1)
  lines(pred,col='red',lwd=2)
}

