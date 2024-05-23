cv.MTPS <- function(xmat, ymat, family, kernel=FALSE, 
                               nfolds = 5,
                               cv = FALSE, residual = TRUE,
                               cv.stacking.nfold = 5, method.step1, method.step2,
                               resid.type=c("deviance", "pearson", "raw"),
                               resid.std=FALSE)
{

  resid.type <- match.arg(resid.type)

  ny <- ncol(ymat)
  # check family input
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

  # metrics
  cindex <- which(family=="gaussian")
  bindex <- which(family=="binomial")
  # if (length(cindex) > 0) {
    metrics.ctn <- matrix(NA, nrow = 1, ncol = length(cindex))
    colnames(metrics.ctn) <- colnames(ymat[, cindex])
    rownames(metrics.ctn) <- "MSE"
    metrics.ctn <- as.data.frame(metrics.ctn)
  # }
  # if (length(bindex) > 0) {
    metrics.bin <- matrix(NA, nrow = 4, ncol = length(bindex))
    colnames(metrics.bin) <- colnames(ymat[, bindex])
    rownames(metrics.bin) <- c("AUC", "Accuracy", "Recall", "precision")
    metrics.bin <- as.data.frame(metrics.bin)
  # }

  idx.cv <- createFolds(rowMeans(xmat), k=nfolds, list=F)
  pred <- ymat; pred[!is.na(pred)] <- NA
  for (i.fold in 1:nfolds) {
    # make train and test data for i.fold-th fold
    y.train <- ymat[idx.cv!=i.fold, ]
    y.test  <- ymat[idx.cv==i.fold, ]
    x.train <- xmat[idx.cv!=i.fold, ]
    x.test  <- xmat[idx.cv==i.fold, ]
    
    fit <- MTPS(xmat = x.train, ymat = y.train, family = family, kernel = kernel,
                           cv = cv, residual = residual,
                           nfold = cv.stacking.nfold,
                           method.step1 = method.step1,
                           method.step2 = method.step2,
                           resid.type = resid.type, resid.std = resid.std)
    pred[idx.cv==i.fold, ] <- predict(fit, x.test)
  }
  # metrics
  if (length(cindex) > 0) {
    metrics.ctn["MSE",] <- apply((pred[, cindex] - ymat[, cindex])^2, 2, mean)
  }
  for (jj in bindex) {
    metrics.bin["AUC", which(jj==bindex)] <- AUC(pred[,jj], outcome=ymat[,jj])
    table <- table((pred[,jj] > 0.5) * 1, ymat[,jj])
    metrics.bin["Accuracy", which(jj==bindex)] <- (table[1,1] + table[2,2]) / sum(table)
    metrics.bin["Recall", which(jj==bindex)] <- table[2,2] / (table[2,2] + table[1,2])
    metrics.bin["precision", which(jj==bindex)] <- table[2,2] / (table[2,2] + table[2,1])
  }
  metrics <- list(continuous = metrics.ctn,
                  binary = metrics.bin)
  return(metrics)
}






