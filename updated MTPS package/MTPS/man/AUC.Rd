\name{AUC}
\alias{AUC}
\title{Area Under Curve}
\usage{
AUC(prob, outcome, cutoff = 1, ROC.plot = FALSE)
}
\arguments{
  \item{prob}{A numeric vector of predicted probability}

  \item{outcome}{A numeric vector of observed binary outcome}

  \item{cutoff}{Number between 0 and 1 to specify where threshold of ROC curve should be truncated. The default value is 1 (no truncation)}

  \item{ROC.plot}{Logical. Whether or not to plot ROC curve}
}
\value{
The value of the area under the curve.
}
\details{

The ROC curve is created by plotting the true positive rate (TPR) against the false positive rate (FPR) at different threshold settings.

By default the total area under the curve is computed, but a truncated AUC statistics can be specified with the cutoff argument. It specifies the bounds of FPR. The common choice of cutoff can be 1 (i.e. no truncate) or 0.2 (i.e. specificity > 0.8)
}
\description{
The AUC function calculates the numeric value of area under the ROC curve (AUC) with the trapezoidal rule and optionally plots the  ROC curve
}

\examples{
set.seed(1)
# simulate predictors
x1 <- rnorm(200)
x2 <- rnorm(200)
# simulate outcome
pr <- 1/(1+exp(-(3 * x1 + 2 * x2 + 1)))
y <- rbinom(200, 1, pr)
df <- data.frame(y = y,x1 = x1, x2 = x2)
# fit logistic regression model on the first 100 observation
lg.model <- glm(y ~ x1 + x2, data = df[1 : 100, ], family="binomial")
# predict outcome for the last 100 observation
prob <- predict(lg.model, df[101:200, c("x1", "x2")], type = "response")
# calculate AUC and plot thr ROC Curve
AUC(prob, y[101:200], ROC=TRUE)
# calculate AUC and plot thr ROC Curve with cutoff
AUC(prob, y[101:200], cutoff=0.2, ROC=TRUE)
}
