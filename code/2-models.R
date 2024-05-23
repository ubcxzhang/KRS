#####models#####
rm(list = ls())
setwd("/Users/caoxiaowen/Desktop/KRS-code-data/")
###load data
load("result/rda files/drug7_gene1067.rda")

library(MTPS)
library(glmnet)
library(RMTL)
library(ggplot2)
library(reshape2)
####prepare x.train/test & kx.train/test####
dat <- c("CCLE","GDSC")
dd <- 2

if(dd==1){
  x.train<-as.matrix(xx1) # if dd==1 indicates CCLE xx1 is training set
  y.train<-as.matrix(log(yy1))
  x.test<-as.matrix(xx2) # GDSC
  y.test<-as.matrix(yy2)
}else{
  x.train<-as.matrix(xx2)
  y.train<-as.matrix(yy2)
  x.test<-as.matrix(xx1)
  y.test<-as.matrix(log(yy1))
}

###normalize-train & test data
for (ii in 1:dim(y.test)[2]) {
  ytestup<-y.test[,ii]-min(y.test[,ii])
  ytestdown<-max(y.test[,ii])-min(y.test[,ii])
  y.test[,ii]<-ytestup/ytestdown
}
for (ii in 1:dim(y.train)[2]) {
  ytestup<-y.train[,ii]-min(y.train[,ii])
  ytestdown<-max(y.train[,ii])-min(y.train[,ii])
  y.train[,ii]<-ytestup/ytestdown
}

### transform x to kernelized x
sigma.train<-mean(as.matrix(dist(x.train)))/sqrt(2)
gamma.train<-1/(2*(sigma.train^2))
kx.train<-exp(-gamma.train*as.matrix(dist(x.train))^2)

k<-matrix(0,dim(x.test)[1],dim(x.train)[1])
for (i in 1:dim(x.test)[1]) { # i:the row number of test data
  for (j in 1:dim(x.train)[1]) { # j:the row number of the train data
    k[i,j]<-sqrt(sum((x.test[i,]-x.train[j,])^2)) # RBF kernel function
    cat(paste0("i:", i,", j: ",j)) # process
  }
}
sigma.test<-mean(k)/sqrt(2)#can change
gamma.test<-1/(2*(sigma.test^2))
kx.test<-exp(-gamma.test*k^2)



####model1 rs####
set.seed(1)
fit.rs <- MTPS(xmat = x.train, ymat = y.train, family = "gaussian",
               cv = FALSE, residual = TRUE,
               method.step1 = glmnet1,
               method.step2 = glmnet.lasso)

pred.rs <- predict(fit.rs, x.test)
mse.rs<-apply((pred.rs-y.test)^2,2, mean)
mse.rs

####model2 krs####
set.seed(1)
fit.krs <- MTPS(xmat = kx.train, ymat = y.train, family = "gaussian",
                cv = FALSE, residual = TRUE,
                method.step1 = glmnet1,
                method.step2 = glmnet.lasso) 
pred.krs <- predict(fit.krs, kx.test)
mse.krs<-apply((pred.krs-y.test)^2,2, mean)
mse.krs

### When users utilize our latest mtps package, they can directly add the option "kernel = TRUE"
# fit.krs <- MTPS(xmat = x.train, ymat = y.train, family = "gaussian",
#                 cv = FALSE, residual = TRUE, kernel = TRUE,
#                 method.step1 = glmnet1,
#                 method.step2 = glmnet.lasso) 
# 
# pred.krs <- predict(fit.krs, x.test)

####model3 l21####
set.seed(1)
y.train.list<-list()
y.test.list<-list()
x.train.list<-list()
x.test.list<-list()

for (i in 1:ncol(y.test)) {
  y.train.list[[i]]<-y.train[,i]
  x.train.list[[i]]<-x.train
  y.test.list[[i]]<-y.test[,i]
  x.test.list[[i]]<-x.test
}

model<-cvMTL(x.train.list, y.train.list, type = "Regression", Regularization = "L21",
             Lam1_seq = seq(0.01, 10, 0.1), Lam2 = 0, G = NULL, k = 2,
             opts = list(init = 0, tol = 10^-3, maxIter = 1000), stratify = FALSE,
             nfolds = 5, ncores = 2, parallel = FALSE)
model1<-MTL(x.train.list, y.train.list, type = "Regression", Regularization = "L21",
            Lam1 = 0.01, Lam1_seq = NULL, Lam2 = 0, opts = list(init = 0, tol
                                                                = 10^-3, maxIter = 1000), G = NULL, k = 2)

pred.L21<-predict(model1,x.test.list)
yhat.L21<-matrix(NA,nrow = dim(y.test)[1],dim(y.test)[2])
for (i in 1:ncol(y.test)) {
  yhat.L21[,i]<-pred.L21[[i]] 
}

mse.L21<-apply((yhat.L21-y.test)^2,2, mean)
mse.L21
median(mse.L21)

####model4 kmtrace####
set.seed(1)
y.train.list<-list()
y.test.list<-list()
x.train.list<-list()
x.test.list<-list()

for (i in 1:ncol(y.test)) {
  y.train.list[[i]]<-y.train[,i]
  x.train.list[[i]]<-x.train
  y.test.list[[i]]<-y.test[,i]
  x.test.list[[i]]<-x.test
}

model<-cvMTL(x.train.list, y.train.list, type = "Regression", Regularization = "Trace",
             Lam1_seq = seq(0.01, 10, 0.1), Lam2 = 0, G = NULL, k = 2,
             opts = list(init = 0, tol = 10^-3, maxIter = 1000), stratify = FALSE,
             nfolds = 5, ncores = 2, parallel = FALSE)


model1<-MTL(x.train.list, y.train.list, type = "Regression", Regularization = "Trace",
            Lam1 = model$Lam1.min, Lam1_seq = NULL, Lam2 = 0, opts = list(init = 0, tol
                                                                          = 10^-3, maxIter = 1000), G = NULL, k = 2)
pred.trace<-predict(model1,x.test.list)
yhat.trace<-matrix(NA,nrow = dim(y.test)[1],dim(y.test)[2])
for (i in 1:ncol(y.test)) {
  yhat.trace[,i]<-pred.trace[[i]] 
}
mse.trace<-apply((yhat.trace-y.test)^2,2, mean)
mse.trace  
median(mse.trace)


####model5 kbmtl####
set.seed(1)

###training function
kbmtl_semisupervised_regression_variational_train <- function(K, Y, parameters) {
  set.seed(parameters$seed)
  
  D <- dim(K)[1]
  N <- dim(K)[2]
  T <- dim(Y)[2]
  R <- parameters$R
  sigma_h <- parameters$sigma_h
  sigma_w <- parameters$sigma_w
  
  Lambda <- list(alpha = matrix(parameters$alpha_lambda + 0.5, D, R), beta = matrix(parameters$beta_lambda, D, R))
  A <- list(mu = matrix(rnorm(D * R), D, R), sigma = array(diag(1, D, D), c(D, D, R)))
  H <- list(mu = matrix(rnorm(R * N), R, N), sigma = array(diag(1, R, R), c(R, R, N)))
  
  epsilon <- list(alpha = matrix(parameters$alpha_epsilon + 0.5 * colSums(!is.na(Y)), T, 1), beta = matrix(parameters$beta_epsilon, T, 1))
  W <- list(mu = matrix(rnorm(R * T), R, T), sigma = array(diag(1, R, R), c(R, R, T)))
  
  KKT <- tcrossprod(K, K)
  
  for (iter in 1:parameters$iteration) {
    # update Lambda
    for (s in 1:R) {
      Lambda$beta[,s] <- 1 / (1 / parameters$beta_lambda + 0.5 * (A$mu[,s]^2 + diag(A$sigma[,,s])))
    }
    # update A
    for (s in 1:R) {
      A$sigma[,,s] <- chol2inv(chol(diag(as.vector(Lambda$alpha[,s] * Lambda$beta[,s]), D, D) + KKT / sigma_h^2))
      A$mu[,s] <- A$sigma[,,s] %*% (tcrossprod(K, H$mu[s,,drop = FALSE]) / sigma_h^2)
    }
    # update H
    for (i in 1:N) {
      indices <- which(is.na(Y[i,]) == FALSE)
      H$sigma[,,i] <- chol2inv(chol(diag(1 / sigma_h^2, R, R) + tcrossprod(W$mu[,indices, drop = FALSE], W$mu[,indices, drop = FALSE] * matrix(epsilon$alpha[indices] * epsilon$beta[indices], R, length(indices), byrow = TRUE)) + apply(W$sigma[,,indices, drop = FALSE] * array(matrix(epsilon$alpha[indices] * epsilon$beta[indices], R * R, length(indices), byrow = TRUE), c(R, R, length(indices))), 1:2, sum)))
      H$mu[,i] <- H$sigma[,,i] %*% (crossprod(A$mu, K[,i]) / sigma_h^2 + tcrossprod(W$mu[,indices, drop = FALSE], Y[i, indices, drop = FALSE] * epsilon$alpha[indices] * epsilon$beta[indices]))
    }
    
    # update epsilon
    for (t in 1:T) {
      indices <- which(is.na(Y[,t]) == FALSE)
      epsilon$beta[t] <- 1 / (1 / parameters$beta_epsilon + 0.5 * (crossprod(Y[indices, t, drop = FALSE], Y[indices, t, drop = FALSE]) - 2 * crossprod(Y[indices, t, drop = FALSE], crossprod(H$mu[,indices, drop = FALSE], W$mu[,t])) + sum((tcrossprod(H$mu[,indices, drop = FALSE], H$mu[,indices, drop = FALSE]) + apply(H$sigma[,,indices, drop = FALSE], 1:2, sum)) * (tcrossprod(W$mu[,t], W$mu[,t]) + W$sigma[,,t]))));
    }
    # update W
    for (t in 1:T) {
      indices <- which(is.na(Y[,t]) == FALSE)
      W$sigma[,,t] <- chol2inv(chol(diag(1 / sigma_w^2, R, R) + epsilon$alpha[t] * epsilon$beta[t] * (tcrossprod(H$mu[,indices, drop = FALSE], H$mu[,indices, drop = FALSE]) + apply(H$sigma[,,indices, drop = FALSE], 1:2, sum))))
      W$mu[,t] <- W$sigma[,,t] %*% (epsilon$alpha[t] * epsilon$beta[t] * H$mu[,indices, drop = FALSE] %*% Y[indices, t, drop = FALSE])
    }
  }
  
  state <- list(Lambda = Lambda, A = A, epsilon = epsilon, W = W, parameters = parameters)
  
}

###predicting function
kbmtl_semisupervised_regression_variational_test <- function(K, state) {
  N <- dim(K)[2]
  T <- dim(state$W$mu)[2]
  
  H <- list(mu = crossprod(state$A$mu, K))
  
  Y <- list(mu = crossprod(H$mu, state$W$mu), sigma = matrix(0, N, T))
  for (t in 1:T) {
    Y$sigma[,t] <- 1 / (state$epsilon$alpha[t] * state$epsilon$beta[t]) + diag(crossprod(H$mu, state$W$sigma[,,t]) %*% H$mu)
  }
  
  prediction <- list(H = H, Y = Y)
}

#initalize the parameters of the algorithm
parameters <- list()

#set the hyperparameters of gamma prior used for projection matrix
parameters$alpha_lambda <- 1
parameters$beta_lambda <- 1

#set the hyperparameters of gamma prior used for output noise
parameters$alpha_epsilon <- 1
parameters$beta_epsilon <- 1

#set the number of iterations
parameters$iteration <- 200

#set the subspace dimensionality
parameters$R <- 20

#set the seed for random number generator used to initalize random variables
parameters$seed <- 1606

#set the standard deviation of hidden representations
parameters$sigma_h <- 0.1

#set the standard deviation of weight parameters
parameters$sigma_w <- 1.0

####predicting section
k<-matrix(0,dim(x.train)[1],dim(x.test)[1])
for (i in 1:dim(x.train)[1]) {
  for (j in 1:dim(x.test)[1]) {
    k[i,j]<-sqrt(sum((x.test[j,]-x.train[i,])^2))
  }
}

sigma.test<-mean(k)/sqrt(2)
gamma.test<-1/(2*(sigma.test^2))
kx.test<-exp(-gamma.test*k^2)

sigma.train<-mean(dist(x.train))/sqrt(2)
gamma.train<-1/(2*(sigma.train^2))
kx.train<-exp(-gamma.train*as.matrix(dist(x.train))^2)

# model
state <- kbmtl_semisupervised_regression_variational_train(kx.train, y.train, parameters)
# predict
pred.kbmtl <- kbmtl_semisupervised_regression_variational_test(kx.test, state)
yhat3<-pred.kbmtl$Y$mu
mse.kbmtl<-apply((yhat3-y.test)^2, 2 ,mean)
mse.kbmtl

####summary mse and absolute bias####
mse5<-cbind(mse.rs,mse.krs,mse.L21,mse.trace,mse.kbmtl)

rsbias <- melt(abs(pred.krs-y.test)-abs(pred.rs-y.test))
l21bias <- melt(abs(pred.krs-y.test)-abs(yhat.L21-y.test))
tracebias <- melt(abs(pred.krs-y.test)-abs(yhat.trace-y.test))
kbmtlbias <- melt(abs(pred.krs-y.test)-abs(yhat3-y.test))

rsbias$Var1 <- rep("RS", length(rsbias$Var1))
l21bias$Var1 <- rep("L21", length(l21bias$Var1))
tracebias$Var1 <- rep("KMTrace", length(tracebias$Var1))
kbmtlbias$Var1 <- rep("KBMTL", length(kbmtlbias$Var1))

t.test(rsbias$value)
t.test(l21bias$value)
t.test(tracebias$value)
t.test(kbmtlbias$value)

bias.abs <- rbind(rsbias, l21bias, tracebias, kbmtlbias)
save(bias.abs,mse5, file = paste0("result/rda files/bias.abs.",dat[dd],"train.rda"))
