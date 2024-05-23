cv.glmnet2 <- function(xx, yy, foldid, alpha=seq(0,10,by=2)/10, lambda=exp(seq(log(10^-8), log(5), length.out=100)),...)
{
  fits <- vector("list",length(alpha))
  names(fits) <- alpha

  for(ii in 1:length(alpha))
  {
    fits[[ii]] <- cv.glmnet(x=xx, y=yy, foldid=foldid, alpha=alpha[ii], lambda=lambda, ...)
  }

  if(length(alpha)==1)
  {
    idx <- 1
  }else
  {
    idx <- which.min(sapply(fits, function(xx) xx$cvm[which(xx$lambda==xx$lambda.1se)]))
  }

  fits[[idx]]$alpha <- alpha[idx]
  return(fits[[idx]])
}
