\name{modify.parameter}
\alias{modify.parameter}
\title{Modify Default Parameters For Base Learner}
\usage{
modify.parameter(FUN, ...)
}
\arguments{
  \item{FUN}{Method}

  \item{...}{Modified arguments}
}
\value{
It returns a new function with modified parameters.
}
\description{
Modify default parameters for methods provided in the package.
}
\examples{
glmnet.lasso <- modify.parameter(glmnet1, alpha=1)
glmnet.ridge <- modify.parameter(glmnet1, alpha=0)
}


