\name{multiFit}
\alias{multiFit}
\title{Fit models on multiple outcomes}
\description{
This function fit individual models to predict each outcome separately.
}
\usage{
multiFit(xmat, ymat, method, family)
}
\arguments{
  \item{xmat}{Matrix of predictors, each row is an observation vector}

  \item{ymat}{Matrix of outcomes. Quantitative for family = "gaussian" and a factor of two levels for family = "binomial"}

  \item{method}{Method for fitting models. It can be one base learner function for all outcomes or a list of base learner functions for each outcome. The list of all base learners can be obtained by \code{list.learners()}.}

  \item{family}{Response type for each response. If all response variable are within the same family it can be "gaussian" or "binomial", otherwise it is a vector of "gaussian" or "binomial" to indicate each response family}
}

\value{
It returns a multiFit object. It is a list of 5 parameters containing information about the fitted models and fitted values for each outcome.
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
fit <- multiFit(xmat = x.train, ymat = y.train,
                method = rpart1, family = "gaussian")
predict(fit, x.test)

# using different base learners for different outcomes
fit.mixOut <- multiFit(xmat = x.train, ymat = y.train,
                method = c(rpart1, rpart1, glmnet.ridge,lm1,lm1),
                family = "gaussian")
predict(fit.mixOut, x.test)
}
