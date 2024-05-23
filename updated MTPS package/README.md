---
title: "An Example of Using KRS of MTPS Package in R."
author: "Xiaowen Cao"
date: "2024-05-11"
output: html_document
---


```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## Introduction

When the relationship between predictors and outcomes is complex and nonlinear, a kernel function can project the data into a higher-dimensional space, simplifying these complex nonlinear relationships. Therefore, we have developed an improved approach based on MTPS, called Kernelized Residual Stacking (KRS), to better capture the nonlinear relationships.

We show an example of using KRS from the article 'Novel machine learning model for predicting cancer drugs' susceptibilities and discovering novel treatments'.

Replication: This repository is designed to guide others in using our tool. If you are interested in the scripts needed to replicate our results, please contact us, and we will provide access to the replication repository. Contact information can be found at the bottom of this page.

## An Example for R Users

This method is currently integrated into the R package, MTPS.

```{r example, eval=FALSE}
if (!require(MTPS)) {
  install.packages("MTPS")
  library(MTPS)
}

set.seed(1)

fit.krs <- MTPS(xmat = x.train, ymat = y.train, family = "gaussian",
                cv = FALSE, residual = TRUE, kernel = TRUE,
                method.step1 = glmnet1,
                method.step2 = glmnet.lasso) 

pred.krs <- predict(fit.krs, x.test)
```

## Contact
[Email me](mailto:xiaowencao@uvic.ca)
