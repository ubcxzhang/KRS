multiFit <- function(xmat, ymat,
                     method,
                     family=family)
{

  ny <- ncol(ymat)
  nx <- ncol(xmat)

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
    stop("each family must be gaussian or binomial")
  }

  # check family method consistency
  if (length(method) == 1) {
    method <- rep(list(method),ny)
  }
  if (length(method) != ny) {
    stop("length of method.step1 must be 1 or the same as response column")
  }
  for (ii in 1:ny) {
    if (!check.match(family[ii], FUN=method[[ii]])) {
      stop("method.step1 must be consistent with response category")
    }
  }

  # initialize matrix to save fitted results
  y.fitted <- ymat; y.fitted[!is.na(y.fitted)] <- NA
  models <- vector("list", ny)
  colnames(y.fitted) <- names(models) <- colnames(ymat)
  fit <- vector("list", ny)

  # colnames of predictor matrix are needed to avoid problem when calling function "predict"
  colnames(xmat) <- paste0("X", 1:nx)

  for (kk in 1:ny)
  {
    fit[[kk]] <- method[[kk]](xmat, ymat[,kk],family=family[kk])
    models[[kk]] <- fit[[kk]]$model
    y.fitted[,kk] <- fit[[kk]]$y.fitted
  }

  multiFit.fits <- list(fit=fit,
              y.fitted=y.fitted,
              model=models,
              method=method,
              family=family)
  class(multiFit.fits) <- "multiFit"
  return(multiFit.fits)
}

predict.multiFit <- function(object, newdata, ...)
{
  ny <- length(object$model)
  pred <- matrix(NA, nrow(newdata), ny)
  mtd <- object$method

  colnames(newdata) <- paste0("X", 1:ncol(newdata))
  for(ii in 1:ny)
  {
    model <- object$model[[ii]]
    pred[,ii] <- object$fit[[ii]]$predFun(model, newdata)

    if(object$family[ii]=="binomial")
    { # predicted are probability should be within [0,1]
      pred[,ii][pred[,ii]>1] <- 1
      pred[,ii][pred[,ii]<0] <- 0
    }
  }
  return(pred)
}



