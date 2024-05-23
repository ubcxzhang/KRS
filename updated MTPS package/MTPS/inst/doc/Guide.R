## ----library, echo = T, results = 'hide', warning=FALSE, message=FALSE--------
library(MTPS)

## ----data---------------------------------------------------------------------
data("HIV")
head(YY)
XX[8:9, 7:10]
dim(YY)
dim(XX)

## ----ctn-split-data, echo = T, results = 'hide'-------------------------------
set.seed(12345)
xmat <- as.matrix(XX)
ymat <- as.matrix(YY)
nobs <- nrow(xmat)
id <- createFolds(rowMeans(XX), k=5, list=F)
training.id <- sample(seq_len(nobs), size = 0.8 * nobs)
y.train <- ymat[training.id, ]
y.test  <- ymat[-training.id, ]
x.train <- xmat[training.id, ]
x.test  <- xmat[-training.id, ]

## ----ctn-noss, echo = T, results = 'hide'-------------------------------------
# no stacking
fit.mult <- multiFit(xmat = x.train, ymat = y.train, method = glmnet.ridge, family = "gaussian")

## ----ctn-train, echo = T, results = 'hide', eval=FALSE------------------------
#  # Standard Stacking
#  fit.ss <- MTPS(xmat = x.train, ymat = y.train, family = "gaussian",
#                              cv = FALSE, residual = FALSE,
#                              method.step1 = glmnet.ridge,
#                              method.step2 = rpart1)
#  # Cross-Validation Stacking
#  fit.cv <- MTPS(xmat = x.train, ymat = y.train, family = "gaussian",
#                              cv = TRUE, residual = FALSE,
#                              method.step1 = glmnet.ridge,
#                              method.step2 = rpart1)
#  # Residual Stacking
#  fit.rs <- MTPS(xmat = x.train, ymat = y.train, family = "gaussian",
#                              cv = FALSE, residual = TRUE,
#                              method.step1 = glmnet.ridge,
#                              method.step2 = rpart1)
#  # Cross-Validation Residual Stacking
#  fit.cvrs <- MTPS(xmat = x.train, ymat = y.train, family = "gaussian",
#                              cv = TRUE, residual = TRUE,
#                              method.step1 = glmnet.ridge,
#                              method.step2 = rpart1)

## ----echo=F,eval=T------------------------------------------------------------
data("Internal")

## ----ctn-predict, echo = T, results = 'hide', eval=FALSE----------------------
#  # no stacking
#  pred.mult <- predict(fit.mult, x.test)
#  # Standard Stacking
#  pred.ss <- predict(fit.ss, x.test)
#  # Cross-Validation Stacking
#  pred.cv <- predict(fit.cv, x.test)
#  # Residual Stacking
#  pred.rs <- predict(fit.rs, x.test)
#  # Cross-Validation Residual Stacking
#  pred.cvrs <- predict(fit.cvrs, x.test)

## ----ctn-outcome, eval=F------------------------------------------------------
#  library(ggplot2)
#  library(reshape2)
#  n.test <- nrow(x.test)
#  ctn.plot.data_matrix <- cbind(rbind(pred.mult, pred.ss, pred.cv, pred.rs, pred.cvrs), y.test[rep(seq_len(n.test), 5), ])
#  ctn.plot.data <-data.frame(rep(c("No-Stacking", "SS", "CV", "RS", "CVRS"), each = n.test),ctn.plot.data_matrix)
#  colnames(ctn.plot.data) <- c("method", paste0(rep("pred.", ncol(y.test)), colnames(y.test)), colnames(y.test))
#  dm1 <- melt(ctn.plot.data[,c("method","ABC","3TC","AZT","D4T", "DDI")], id=c("method"))
#  dm2 <- melt(ctn.plot.data[,c("method","pred.ABC","pred.3TC","pred.AZT","pred.D4T", "pred.DDI")], id=c("method"))
#  ctn.plot.data <- cbind(dm1, dm2[, -1])
#  colnames(ctn.plot.data) <- c("method", "Y", "yVal", "predict", "predictVal")
#  ctn.plot.data$method <- factor(ctn.plot.data$method, unique(as.character(ctn.plot.data$method)))
#  ctn.plot.data$yVal <- as.numeric(ctn.plot.data$yVal)
#  ctn.plot.data$predictVal <- as.numeric(ctn.plot.data$predictVal)
#  ggplot(ctn.plot.data) +
#    geom_point(aes(x=predictVal, y=yVal, color=method), size = 0.5, alpha = 0.8) +
#    geom_abline(slope=1, alpha = 0.2) +
#    coord_fixed() +
#    ylab("Testing Data Outcome") + xlab("Predicted Outcome on Testing Data") +
#    scale_x_discrete(breaks = NULL) +
#    scale_y_discrete(breaks = NULL) +
#    theme_bw() +
#    theme(axis.text = element_blank(),
#          strip.placement = "outside",
#          strip.background = element_blank()) +
#    facet_grid(Y~method)

## ----warning=FALSE, echo=FALSE------------------------------------------------
library(ggplot2)
library(reshape2)
ggplot(ctn.plot.data) +
  geom_point(aes(x=predictVal, y=yVal, color=method), size = 0.5, alpha = 0.8) +
  geom_abline(slope=1, alpha = 0.2) +
  coord_fixed() +
  ylab("Testing Data Outcome") + xlab("Predicted Outcome on Testing Data") +
  scale_x_discrete(breaks = NULL) +
  scale_y_discrete(breaks = NULL) +
  theme_bw() +
  theme(axis.text = element_blank(),
        strip.placement = "outside",
        strip.background = element_blank()) +
  facet_grid(Y~method)

## ----bin-data, echo = T, results = 'hide', eval=FALSE-------------------------
#  # https://hivdb.stanford.edu/pages/published_analysis/genophenoPNAS2006/CUTOFFS/drug.cutoffs
#  # cutoff value to be used to define drug resistent
#  cutoffs <- c(2,3,3,1.5,1.5)
#  ymat.bin <- ymat
#  xmat.bin <- xmat
#  for(ii in 1:5) ymat.bin[,ii] <- (10^ymat[,ii] < cutoffs[ii]) * 1
#  y.train.bin <- ymat.bin[training.id, ]
#  y.test.bin  <- ymat.bin[-training.id, ]
#  x.train.bin <- xmat.bin[training.id, ]
#  x.test.bin  <- xmat.bin[-training.id, ]

## ----bin-train, echo = T, results = 'hide', eval=FALSE------------------------
#  fit.prs.std <- MTPS(xmat = x.train.bin, ymat=y.train.bin,
#                                 family = "binomial",
#                                 cv = FALSE, residual = TRUE,
#                                 method.step1 = rpart1,
#                                 method.step2 = lm1,
#                                 resid.type = "pearson", resid.std = TRUE)
#  pred.prs.std <- predict(fit.prs.std, x.test.bin)

## ----bin-outcome, echo = T----------------------------------------------------
for (yy in 1 : ncol(y.test.bin)) {
  print(colnames(y.test.bin)[yy])
  print(table((pred.prs.std[,yy] > 0.5) * 1, y.test.bin[,yy]))
}


## ----mix-data, echo = T, results = 'hide', eval=FALSE-------------------------
#  ymat.mix <- cbind(ymat[,1:3], ymat.bin[,4:5])
#  xmat.mix <- xmat
#  y.train.mix <- ymat.mix[training.id, ]
#  y.test.mix  <- ymat.mix[-training.id, ]
#  x.train.mix <- xmat.mix[training.id, ]
#  x.test.mix  <- xmat.mix[-training.id, ]

## ----mix-training, echo = T, results = 'hide', warning=FALSE, message=FALSE, eval=FALSE----
#  fit.mix.rs <- MTPS(xmat = x.train.mix, ymat = y.train.mix,
#                  family=c("gaussian","gaussian","gaussian","binomial","binomial"),
#                  cv = FALSE, residual = TRUE,
#                  method.step1 = glmnet.lasso,
#                  method.step2 = rpart1)
#  pred.mix.rs <- predict(fit.mix.rs, x.test.mix)

## ----mix-outcome, eval=FALSE--------------------------------------------------
#  n.test <- nrow(x.test)
#  mix.plot.data <- cbind(rep(colnames(y.test.mix)[1:3], each=nrow(y.test.mix)),
#                         rbind(cbind(pred.mix.rs[, 1], y.test.mix[, 1]),
#                               cbind(pred.mix.rs[, 2], y.test.mix[, 2]),
#                               cbind(pred.mix.rs[, 3], y.test.mix[, 3])))
#  colnames(mix.plot.data) <- c("Y", "predict", "outcome")
#  mix.plot.data <- as.data.frame(mix.plot.data)
#  mix.plot.data$predict <- as.numeric(as.character(mix.plot.data$predict))
#  mix.plot.data$outcome <- as.numeric(as.character(mix.plot.data$outcome))
#  ggplot(mix.plot.data) +
#    geom_point(aes(x=predict, y=outcome, color=Y), size = 0.5, alpha = 0.8) +
#    ylab("Outcome of Testing Data") + xlab("Predicted Outcome of Testing Data") +
#    scale_x_discrete(breaks = NULL) +
#    scale_y_discrete(breaks = NULL) +
#    geom_abline(slope=1, alpha = 0.2) +
#    coord_fixed() +
#    theme_bw() +
#    theme(legend.title = element_blank(),
#          axis.text = element_blank(),
#          strip.placement = "outside",
#          strip.background = element_blank()) +
#    facet_grid(~Y)

## ----echo=FALSE---------------------------------------------------------------
ggplot(mix.plot.data) +
  geom_point(aes(x=predict, y=outcome, color=Y), size = 0.5, alpha = 0.8) +
  ylab("Outcome of Testing Data") + xlab("Predicted Outcome of Testing Data") +
  scale_x_discrete(breaks = NULL) +
  scale_y_discrete(breaks = NULL) +
  geom_abline(slope=1, alpha = 0.2) +
  coord_fixed() +
  theme_bw() +
  theme(legend.title = element_blank(),
        axis.text = element_blank(),
        strip.placement = "outside",
        strip.background = element_blank()) +
  facet_grid(~Y)

## ----echo=T, eval=T-----------------------------------------------------------
for (yy in 4 : 5) {
  print(colnames(y.test.mix)[yy])
  print(table((pred.mix.rs[,yy] > 0.5) * 1, y.test.mix[,yy]))
}

## ----mix-mtd, echo = T, warning=FALSE,message=FALSE, eval=FALSE---------------
#  fit.mixOut <- MTPS(xmat=x.train, ymat=y.train, family="gaussian",
#                  method.step1 =
#                    c(glmnet.lasso,glmnet.lasso,glmnet.lasso,lm1,lm1),
#                  method.step2 =
#                    c(rpart1,glmnet.lasso,glmnet.lasso,glmnet.lasso,glmnet.lasso))
#  pred <- predict(fit.mixOut, x.test)

## ----method-mod, eval=F, echo=T,eval=FALSE------------------------------------
#  glmnet.lasso <- modify.parameter (glmnet1, alpha=1)
#  glmnet.ridge <- modify.parameter (glmnet1, alpha=0)

## ----method-new, eval=F, echo=T,eval=FALSE------------------------------------
#  glm1 <- function(xmat, ymat, family, ...) {
#    tmp0 <- data.frame(yy=ymat, xmat)
#    model <- glm(yy~., data=tmp0, family=family, ...)
#    y.fitted <- fitted(model)
#    predFun <- function(model,xnew){
#      predict(model, newdata=data.frame(xnew), type="response")
#    }
#    return(list(model=model,y.fitted=y.fitted, predFun=predFun))
#  }

## ----learner-compare, echo=T,eval=F-------------------------------------------
#  nsim <- 20
#  mse.lasso.lm <- matrix(NA, nrow = nsim, ncol = ncol(ymat))
#  mse.lasso.lm <- as.data.frame(mse.lasso.lm)
#  colnames(mse.lasso.lm) <-colnames(ymat)
#  mse.ridge.lm <- mse.lasso.lm
#  for (ii in 1:nsim) {
#    set.seed(ii)
#    # lasso stacking with lm
#    mse.lasso.lm[ii,] <- cv.MTPS(xmat, ymat, family="gaussian",
#                              cv = FALSE, residual = TRUE,
#                              method.step1=glmnet.lasso, method.step2=lm1,
#                              resid.std=FALSE)$continuous
#    # ridge stacking with lm
#    mse.ridge.lm[ii,] <- cv.MTPS(xmat, ymat, family="gaussian",
#                              cv = FALSE, residual = TRUE,
#                              method.step1=glmnet.ridge, method.step2=lm1,
#                              resid.std=FALSE)$continuous
#  }
#  # box plot
#  mse.data <- data.frame(lasso=apply(mse.lasso.lm,1,mean),
#                         ridge=apply(mse.ridge.lm,1,mean))
#  mse.data <- melt(mse.data,measure.vars = c("lasso","ridge"))
#  
#  colnames(mse.data) <- c("Learner", "MSE")
#  ggplot(mse.data) +
#    geom_boxplot(aes(x=Learner, y=MSE, fill=Learner)) +
#    ggtitle("Boxplot of average Mean Square Error") +
#    theme_bw() +
#    theme(legend.position = "none",
#          plot.title = element_text(hjust = 0.5),
#          axis.title.x = element_blank())

## ----echo=F,eval=T------------------------------------------------------------
# box plot
ggplot(mse.data) +
  geom_boxplot(aes(x=Learner, y=MSE, fill=Learner)) +
  ggtitle("Boxplot of average Mean Square Error") +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank())

