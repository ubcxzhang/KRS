\name{MTPS-internal}
\alias{MTPS-internal}
\alias{trapezoid}
\alias{check.match}
\alias{createFolds}
\alias{cv.glmnet2}
\alias{cv.multiFit}
\alias{max_finite}
\alias{min_finite}
\alias{resid.bin}
\alias{rs.multiFit}
\alias{predict.cv.multiFit}
\alias{predict.rs.multiFit}
\title{Internal functions}
\description{
Internal functions for MTPS package
}
\usage{
trapezoid(FPR,TPR,cc)

check.match(family, FUN)

createFolds(y, k = 10, list = TRUE, returnTrain = FALSE)

cv.glmnet2(xx, yy, foldid, alpha=seq(0,10,by=2)/10,
lambda=exp(seq(log(10^-8), log(5), length.out=100)),...)

cv.multiFit(xmat, ymat, nfold=5, method, family=family)

max_finite(xx)

min_finite(xx)

resid.bin(ymat, yhat, xmat=NULL, type=c("deviance", "pearson", "raw"), resid.std=F)

rs.multiFit(yhat, ymat, xmat=NULL, family,
            resid.type=c("deviance", "pearson", "raw"), resid.std=F,
            method)
}
\details{These are not intended for use by users.
\code{check.match} check whether a given method matches with the response family.
\code{createFolds} the createFolds function is cited from the 'caret' library, which is a large r package for machine learning. To improve the independency of our package we copied this single function instead of loading the whole caret package. It is safely to ignore the warning if you have loaded the caret the package.
\code{cv.multiFit} used to fit models in cross-validation stacking
\code{rs.multiFit} used to fit models in residual stacking
\code{resid.bin} calculate residual of different types
}
\keyword{internal}
