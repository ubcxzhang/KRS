\name{MTPS}
\alias{MTPS}
\title{Fit Models using Revised Stacking Algorithm}
\usage{
MTPS(xmat, ymat, family, kernel = FALSE,
  cv = FALSE, residual = TRUE, nfold = 5,
  method.step1, method.step2,
  resid.type = c("deviance", "pearson", "raw"), resid.std = FALSE)
}
\arguments{
  \item{xmat}{Predictor matrix, each row is an observation vector}

  \item{ymat}{Responses matrix. Quantitative for family = "gaussian" and a factor of two levels for family = "binomial"}

  \item{family}{Response type for each response. If all response variable are within the same family it can be "gaussian" or "binomial", otherwise it is a vector with elements "gaussian" and "binomial" to indicate each response family}

  \item{kernel}{Logical, indicate if use RBF kernel functions on predictor matrix}
  
  \item{cv}{Logical, indicate if use Cross-Validation Stacking algorithm}

  \item{residual}{Logical, indicate if use Residual Stacking algorithm}

  \item{nfold}{Integer, number of folds for Cross-Validation Stacking algorithm. The default value is 5}

  \item{method.step1}{Base Learners for fitting models in Step 1 of Stacking Algorithm. It can be one base learner function for all outcomes or a list of base learner functions for each outcome. The list of all base learners can be obtained by \code{list.learners()}}

  \item{method.step2}{Base Learners for fitting models in Step 2 of Stacking Algorithm. (see above)}

  \item{resid.type}{The residual type for Residual Stacking}

  \item{resid.std}{Logical, whether or not use standardized residual}

}
\value{
It returns a MTPS object. It is a list of 4 parameters containing information about step 1 and step 2 models and the revised stacking algorithm method.
}
\description{
Fit a model using standard stacking algorithm or revised stacking algorithms to simultaneous predicte multiple outcomes
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

# Residual Stacking
fit.rs <- MTPS(xmat = x.train, ymat = y.train,
  family = "gaussian",cv = FALSE, residual = TRUE,
  method.step1 = rpart1, method.step2 = lm1)
predict(fit.rs, x.test)

# using different base learners for different outcomes
 fit.mixOut <- MTPS(xmat=x.train, ymat=y.train,
  family="gaussian",cv = FALSE, residual = TRUE,
  method.step1 = c(rpart1,glmnet.ridge,rpart1,lm1,lm1),
  method.step2 = c(rpart1,lm1,lm1,lm1, glmnet.ridge))
predict(fit.mixOut, x.test)
}


