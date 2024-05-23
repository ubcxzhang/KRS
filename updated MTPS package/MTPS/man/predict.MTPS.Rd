\name{predict.MTPS}
\alias{predict.MTPS}
\title{Make predictions from a "MTPS" model}
\description{
This function makes predictions from a revised stacking model.}
\usage{
\method{predict}{MTPS}(object, newdata, ...)
}
\arguments{
  \item{object}{A fitted object from \code{"MTPS"}}

  \item{newdata}{Matrix of new predictors at which predictions are to be made}
  
  \item{...}{additional arguments affecting the predictions produced}
}
\value{
  The predicted value from new predictors.
}

\examples{
data("HIV")
set.seed(1)
xmat <- as.matrix(XX)
ymat <- as.matrix(YY)
id <- createFolds(rowMeans(XX), k=5, list=FALSE)
training.id <- id != 1
y.train <- ymat[training.id, ]
y.test  <- ymat[!training.id, ]
x.train <- xmat[training.id, ]
x.test  <- xmat[!training.id, ]
# Cross-Validation Residual Stacking
fit.rs <- MTPS(xmat = x.train, ymat = y.train,
  family = "gaussian",cv = FALSE, residual = TRUE,
  method.step1 = rpart1, method.step2 = lm1)
pred.rs <- predict(fit.rs, x.test)
}
