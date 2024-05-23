check.match <- function(family, FUN) {

  continuous <- c(lm1,glm1,glmnet1,rpart1,glmnet.lasso,glmnet.ridge)
  binary <- c(glm1,glmnet1,rpart1,glmnet.lasso,glmnet.ridge,lda1,qda1,knn1,svm1)

  clength <- length(continuous)
  blength <- length(binary)

  userDefine <- TRUE

  for (ii in 1:clength) {
    if (identical(FUN, continuous[[ii]]) ) {
      userDefine <- FALSE
      if (family == "gaussian") {
        return(TRUE)
      }
    }
  }
  for (ii in 1:blength) {
    if (identical(FUN, binary[[ii]]) ) {
      userDefine <- FALSE
      if (family == "binomial") {
        return(TRUE)
      }
    }
  }

# 	if (userDefine){
# 		warning("Please make sure that method is consistnet with outcome type")
# 		return(TRUE)
# 	}

  return(FALSE)

}

