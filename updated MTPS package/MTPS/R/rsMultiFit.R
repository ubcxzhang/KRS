rs.multiFit <- function(yhat, ymat, xmat = NULL, kernel = NULL,
                        family,
                        resid.type=c("deviance", "pearson", "raw"), resid.std = F,
                        method)
{

  resid.type <- match.arg(resid.type)

  ny <- ncol(ymat)

  if (length(family) == 1) {
    if (!family %in% c("gaussian", "binomial")) {
      stop("family must be gaussian or binomial")
    }
    if (family == "gaussian") {
      family = rep("gaussian", ny)
    } else if (family == "binomial") {
      family = rep("binomial", ny)
    }
  }

  if (length(family) != ny) {
    stop("length of family must be consistent with response")
  }
  if (sum(family %in% c("gaussian", "binomial")) != ny) {
    stop("family must be gaussian or binomial or their combination")
  }

  if (length(method) == 1) {
    method <- rep(list(method),ny)
  }
  if (length(method) != ny) {
    stop("length of method must be 1 or the same as response column")
  }
  for (ii in 1:ny) {
    if (!check.match("gaussian", FUN=method[[ii]])) {
      stop("residual stacking does not allow binary outcome model in second step")
    }
  }

  y.fitted <- ymat; y.fitted[!is.na(y.fitted)] <- NA
  models <- vector("list", ny)
  colnames(y.fitted) <- names(models) <- colnames(ymat)
  fit <- vector("list", ny)

  # colnames of predictor matrix are needed to avoid problem when calling function "predict"
  colnames(yhat) <- paste0("X", 1:ny)

  # derive residuals to be predicted by other outcome variables

  for (kk in 1:ny)
  {
    if(family[kk]=="gaussian") resi <- ymat-yhat
    if(family[kk]== "binomial") resi <- resid.bin(ymat, yhat, xmat, type=resid.type, resid.std=resid.std)
    xmat=yhat[,-kk]
    if(ny==2){
      xmat <-as.data.frame(xmat)
      names(xmat)<-paste0("X",c(1:ny)[-kk])
    }

    fit[[kk]] <- method[[kk]](xmat=xmat, ymat=resi[,kk], family="gaussian")
    models[[kk]] <- fit[[kk]]$model
    y.fitted[,kk] <- fit[[kk]]$y.fitted
  }

  rsMultiFit.fits <- list(fit=fit,
                     y.fitted=y.fitted,
                     family=family,
                     models=models,
                     method=method,
                     resid.type=resid.type,
                     resid.std=resid.std)
  class(rsMultiFit.fits) <- "rs.multiFit"
  return(rsMultiFit.fits)
}

predict.rs.multiFit <- function(object, newdata, ...) {

  ny <- ncol(object$y.fitted)
  resid.mat <- newdata; resid.mat[!is.na(resid.mat)] <- NA
  mtd <- object$method

  pred<- newdata; pred[!is.na(pred)] <- NA

  colnames(newdata) <- paste0("X", 1:ny)
  for(ii in 1:ny)
  {
    xx <- newdata[,-ii]
    if(ny==2){
      xx <-as.data.frame(xx)
      names(xx)<-paste0("X",c(1:ny)[-ii])
    }
    
    model <- object$model[[ii]]
    resid.mat[,ii] <- object$fit[[ii]]$predFun(model, xx)

  }

  cindex<-object$family=="gaussian"
  bindex<-object$family=="binomial"

  # continuous outcome
  pred[,cindex] <- newdata[,cindex] + resid.mat[,cindex]
  # binary outcome
  if (sum(bindex) != 0) {
    if(object$resid.std)
    {
      QQ <- qr.Q(qr(newdata[,bindex]))  # Q matrix of QR decomposition
      hats <- rowSums(QQ^2) # diagonal of hat matrix X(X'X)^(-1)X', used to standardize residuals
    }else
    {
      hats <- rep(0, nrow(matrix(newdata[,bindex])))
    }
    if(object$resid.type == "pearson")
      pred[,bindex] <- newdata[,bindex] + resid.mat[,bindex]*sqrt(newdata[,bindex]*(1-newdata[,bindex]))*sqrt(1-hats)
    if(object$resid.type == "raw")
      pred[,bindex] <- newdata[,bindex] + resid.mat[,bindex]
    if(object$resid.type == "deviance")
    {
      dev0 <- -sqrt(-2*log(1-newdata[,bindex]))
      dev1 <-  sqrt(-2*log(newdata[,bindex]))
      # avoid infinity residuals
      m1 <- max_finite(dev1)
      m2 <- min_finite(dev0)
      dev1[dev1 > m1] <- m1+100
      dev0[dev0 < m2] <- m2-100

      res0 <- abs(resid.mat[,bindex]-dev0)
      res1 <- abs(resid.mat[,bindex]-dev1)
      pred[,bindex] <- res0/(res0+res1)
    }

    pred[,bindex][pred[,bindex]>1] <- 1
    pred[,bindex][pred[,bindex]<0] <- 0
  }

  return(pred)
}


