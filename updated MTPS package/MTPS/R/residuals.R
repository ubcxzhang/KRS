max_finite <- function(xx) max(xx[is.finite(xx)])
min_finite <- function(xx) min(xx[is.finite(xx)])

# calculate residuals for binary predictions, based on observed binary values and predicted probabilities
resid.bin <- function(ymat, yhat, xmat=NULL,
                      type=c("deviance", "pearson", "raw"),
                      resid.std=F)
{
  if(resid.std)  # obtain the leverage h_k used to standardize residuals
  {
    if(is.null(xmat)) stop("standardized residual need info about hat matrix")
    QQ <- qr.Q(qr(xmat))  # Q matrix of QR decomposition
    hats <- rowSums(QQ^2) # diagonal of hat matrix X(X'X)^(-1)X', used to standardize residuals
  }else
  {
    hats <- 0
  }

  resi0 <- ymat-yhat # raw residual (for continuous outcome models)
  if(type=="deviance") resi <- ifelse(ymat==1, sqrt(-2*log(yhat)), -sqrt(-2*(log(1 - yhat))))
  if(type=="pearson")  resi <- ifelse(ymat==1, exp(log(1-yhat)/2-log(yhat)/2), -exp(log(yhat)/2-log(1-yhat)/2))
  if(type=="raw") 		 resi <- ymat-yhat
  if(resid.std ) 			 resi <- resi/sqrt(1-hats)

  # avoid infinity residuals
  m1 <- max_finite(resi)
  m2 <- min_finite(resi)
  resi[resi > m1] <- m1+100
  resi[resi < m2] <- m2-100
  return(resi)
}
