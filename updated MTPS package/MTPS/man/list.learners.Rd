\name{list.learners}
\alias{list.learners}
\alias{KNN1}
\alias{glm1}
\alias{glmnet.lasso}
\alias{glmnet.ridge}
\alias{glmnet1}
\alias{lda1}
\alias{lm1}
\alias{qda1}
\alias{rpart1}
\alias{svm1}
\title{List Available Base Learners}
\description{
This function lists all base learners provided in the package.
}

\usage{
list.learners()
}
\value{
The name of all base learners provided in the package
}
\details{
lm1: linear regression

glm1: generalized linear models

glmnet1: Does k-fold cross-validation to chose best alpha and lambda for generalized linear models via penalized maximum likelihood.

glmnet.lasso: LASSO, lambda is chose by k-fold cross-validation for glmnet

glmnet.ridge: Ridge regression, lambda is chose by k-fold cross-validation for glmnet

rpart1: regression tree

lda1: linear discriminant analysis

qda1: quadratic discriminant analysis

KNN1: k-nearest neighbour classification, k is chose by cross-validation

svm1: support vector machine

}

\examples{
list.learners()
}
