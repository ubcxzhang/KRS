\name{cv.MTPS}
\alias{cv.MTPS}
\title{Evaluation using Cross-Validation}
\usage{
cv.MTPS(xmat, ymat, family, nfolds = 5, kernel= FALSE,
                   cv = FALSE, residual = TRUE,
                   cv.stacking.nfold = 5, method.step1, method.step2,
                   resid.type=c("deviance", "pearson", "raw"),
                   resid.std=FALSE)

}
\arguments{
  \item{xmat}{Predictor matrix, each row is an observation vector}

  \item{ymat}{Responses matrix. Quantitative for family = "gaussian" and a factor of two levels for family = "binomial"}

  \item{family}{Response type for each response. If all response variable are within the same family it can be "gaussian" or "binomial", otherwise it is a vector with elements "gaussian" and "binomial" to indicate each response family}
  
 \item{kernel}{Logical, indicate if use RBF kernel functions on predictor matrix}
 
  \item{nfolds}{Integer, number of folds for Cross-Validation to evaluate the performance of stacking algorithms.}

  \item{cv}{Logical, indicate if use Cross-Validation Stacking algorithm}

  \item{residual}{Logical, indicate if use Residual Stacking algorithm}

  \item{cv.stacking.nfold}{Integer, number of folds for Cross-Validation Stacking algorithm. The default value is 5}

  \item{method.step1}{Base Learners for fitting models in Step 1 of Stacking Algorithm. It can be one base learner function for all outcomes or a list of base learner functions for each outcome. The list of all base learners can be obtained by \code{list.learners()}}

  \item{method.step2}{Base Learners for fitting models in Step 2 of Stacking Algorithm. (see above)}

  \item{resid.type}{The residual type for Residual Stacking}

  \item{resid.std}{Logical, whether or not use standardized residual}

}
\value{
It returns the mean squared error of continuous outcomes. AUC, accuracy, recall and precision for binary outcomes of predictions using cross-validation.
}
\description{
Use cross-validation to evaluate model performance.
}
\examples{
data("HIV")
cv.MTPS(xmat=XX, ymat=YY, family="gaussian", nfolds=2,
        method.step1=rpart1, method.step2=lm1)
}

