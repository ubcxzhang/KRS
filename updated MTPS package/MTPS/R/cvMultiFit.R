cv.multiFit <- function(xmat, ymat,
                        nfold=5,
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

  ymat.stand <- scale(ymat)
  km1 <- kmeans(ymat.stand, 10, nstart=100)
  # km1 <- kmeans(ymat, 10, nstart=100)
  idx.cv <- createFolds(factor(km1$cluster), k=nfold, list=F)

  # prepare results to be returned
  y.fitted <- ymat; y.fitted[!is.na(y.fitted)] <- NA
  models <- vector("list", nfold)
  names(models) <- paste0("fold",1:nfold)
  fits <- vector("list", nfold)

  for(ii in 1:nfold){
    #make train and test data for ii-th fold
    y.train <- ymat[idx.cv!=ii, ]
    y.test <-  ymat[idx.cv==ii, ]
    x.train <- xmat[idx.cv!=ii, ]
    x.test <-  xmat[idx.cv==ii, ]
    colnames(x.test) <- paste0("X",1:ncol(x.test))

    ###################### ??? should use train data
    # fits[[ii]] <- multiFit(xmat=xmat, ymat=ymat,method,family)

    fits[[ii]] <- multiFit(xmat=x.train, ymat=y.train,method,family)

    y.fitted[idx.cv==ii,] <- predict.multiFit(fits[[ii]], x.test)
    models[[ii]] <- fits[[ii]]$model
  }

  cv.multiFit.fits <- list(fit=fits,
                        y.fitted=y.fitted,
                        model=models,
                        method=method,
                        family=family)
  class(cv.multiFit.fits) <- "cv.multiFit"
  return(cv.multiFit.fits)
  # fit <- list(fit=fits,
  #             y.fitted=y.fitted,
  #             model=models,
  #             method=method,
  #             family=family)
  #
  # return(fit)
}



predict.cv.multiFit <- function(object, newdata, ...)
{
  ny <- ncol(object$y.fitted)
  nfold <- length(object$model)
  temp0 <- array(NA, c(nrow(newdata), ny, nfold))
  clas <- object$method

  colnames(newdata) <- paste0("X", 1:ncol(newdata))
  for(ii in 1:ny) for(jj in 1:nfold)
  {
    model <- object$model[[jj]][[ii]]
    temp0[,ii,jj] <- object$fit[[jj]][["fit"]][[ii]]$predFun(model, newdata)
  }
  pred <- apply(temp0, c(1,2), mean)

  bindex<-object$family=="binomial"
  # predicted are probability should be within [0,1] for binary outcome
  pred[,bindex][pred[,bindex]>1] <- 1
  pred[,bindex][pred[,bindex]<0] <- 0

  return(pred)
}






