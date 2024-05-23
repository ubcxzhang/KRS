modify.parameter <- function(FUN, ...) {
  if (!is.function(FUN)) {
    stop("FUN must be a valid function")
  }
  .FUN <- FUN
  args <- list(...)
  invisible(lapply(seq_along(args), function(i) {
    formals(.FUN)[[names(args)[i]]] <<- args[[i]]
  }))
  .FUN
}

lm1 <- function(xmat, ymat, family, ...) {
  tmp0 <- data.frame(yy=ymat, xmat)
  model <- lm(yy~., data=tmp0, ...)
  y.fitted <- fitted(model)
  predFun <- function(model,xnew){
    predict(model, newdata=data.frame(xnew), type="response")
  }
  return(list(model=model,y.fitted=y.fitted, predFun=predFun))
}

glm1 <- function(xmat, ymat, family, ...) {
  if(family=="binomial") ymat <- factor(ymat)
  tmp0 <- data.frame(yy=ymat, xmat)
  model <- glm(yy~., data=tmp0, family=family, ...)
  y.fitted <- fitted(model)
  predFun <- function(model,xnew){
    predict(model, newdata=data.frame(xnew), type="response")
  }
  return(list(model=model,y.fitted=y.fitted, predFun=predFun))
}

glmnet1 <- function(xmat, ymat,family,alpha=seq(0,10,by=2)/10, ...) {
  tmp0 <- data.frame(yy=ymat, xmat)
  if(family=="binomial") ymat <- factor(ymat)
  foldid <- createFolds(ymat, k=5, list=F)
  model <-  cv.glmnet2(xmat, ymat, alpha=alpha, foldid=foldid, family=family, ...)
  coef.mat <- as.numeric(coef(model,s="lambda.1se"))
  y.fitted <- predict(model, xmat, s="lambda.1se", type="response")
  predFun <- function(model,xnew){
    predict(model, newx=as.matrix(xnew), s="lambda.1se", type="response")
  }
  return(list(model=model,y.fitted=y.fitted, predFun=predFun))
}

glmnet.lasso <- modify.parameter (glmnet1, alpha=1)
glmnet.ridge <- modify.parameter (glmnet1, alpha=0)

rpart1 <- function(xmat, ymat, family, ...){
  tmp0 <- data.frame(yy=ymat, xmat)
  fit0 <- rpart(yy~., data=tmp0, ...)
  model <- prune(fit0, cp=fit0$cptable[which.min(fit0$cptable[,"xerror"]),"CP"])
  y.fitted <- predict(model, newdata=data.frame(xmat))
  predFun <- function(model,xnew){
    predict(model, newdata=data.frame(xnew))
  }
  return(list(model=model,y.fitted=y.fitted, predFun=predFun))
}

lda1 <- function(xmat, ymat, family, ...){
  if(family=="binomial") ymat <- factor(ymat)
  tmp0 <- data.frame(yy=ymat, xmat)
  model <- lda(yy~., data=tmp0, ...)
  y.fitted <- predict(model, newdata=data.frame(xmat))$posterio[,"1"]
  predFun <- function(model,xnew){
    predict(model, newdata=data.frame(xnew))$posterio[,"1"]
  }
  return(list(model=model,y.fitted=y.fitted, predFun=predFun))
}

qda1 <- function(xmat, ymat, family, ...){
  if(family=="binomial") ymat <- factor(ymat)
  tmp0 <- data.frame(yy=ymat, xmat)
  model <- qda(yy~., data=tmp0, ...)
  y.fitted <- predict(model, newdata=data.frame(xmat))$posterio[,"1"]
  predFun <- function(model,xnew){
    predict(model, newdata=data.frame(xnew))$posterio[,"1"]
  }
  return(list(model=model,y.fitted=y.fitted, predFun=predFun))
}

KNN1 <- function(xmat, ymat, family, ...){
  if(family=="binomial") ymat <- factor(ymat)
  tmp0 <- data.frame(yy=ymat, xmat)
  # Select k applying cross-validation
  knn_acc <- rep(NA, sqrt(length(ymat)))
  for (i in 1:sqrt(length(ymat))){
    # print(i)
    results<-knn.cv(train=xmat, cl=ymat, k = i, prob = TRUE)
    knn_acc[i]<-(table(results,ymat)[1,1]+table(results,ymat)[2,2])/(sum(table(results,ymat)))
  }
  k <- which.max(knn_acc)
  model <- list(model=knn.cv(train=xmat, cl=ymat, k = k, prob = TRUE),
                  train=xmat,cl=ymat,k=k)
  y.fitted <- as.numeric(as.character(model$model)) * attributes(model$model)$prob +
    (1-as.numeric(as.character(model$model))) * (1-attributes(model$model)$prob)
  predFun <- function(model,xnew){
    knn.fit <- knn(train=model$train, test=xnew, cl=model$cl, k=model$k, prob = TRUE)
    as.numeric(as.character(knn.fit)) * attributes(knn.fit)$prob +
      (1-as.numeric(as.character(knn.fit))) * (1-attributes(knn.fit)$prob)
  }
  return(list(model=model,y.fitted=y.fitted, predFun=predFun))
}

svm1 <- function(xmat, ymat, family, kernel = "linear", ...){
  if(family=="binomial") ymat <- factor(ymat)
  tmp0 <- data.frame(yy=ymat, xmat)
  model <- svm(yy~.,data=tmp0, kernel = kernel, probability = TRUE, ...)
  pred <- predict(model, xmat, probability=TRUE)
  y.fitted <- attr(pred, "probabilities")[, "1"]
  predFun <- function(model,xnew){
    pred <- predict(model, xnew, probability=TRUE)
    attr(pred, "probabilities")[, "1"]
  }
  return(list(model=model,y.fitted=y.fitted, predFun=predFun))
}
