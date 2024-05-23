list.learners <- function(){

  learners <- c("lm1", "glm1", "glmnet1", "glmnet.lasso", "glmnet.ridge", "rpart1",
                "lda1", "qda1", "KNN1", "svm1")
  print("Models that can be chosen are:")
  print(learners)

}
