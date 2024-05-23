trapezoid<- function(FPR,TPR,cc)
{
	ordFPR=order(FPR)
	FPR=FPR[ordFPR]
	TPR=TPR[ordFPR]
	if (max(FPR)<cc)
	{
		FPRcut=FPR
		TPRcut=TPR
	}else if(min(FPR)>cc)
	{
		stop("threshold smaller than smallest allowed value")
	}else
	{
		FPRlcc=sum(FPR<cc)
		FPR1=FPR[FPRlcc]; FPR2=FPR[FPRlcc+1]
		TPR1=FPR[FPRlcc]; TPR2=TPR[FPRlcc+1]
		TPRcc=( TPR1*(cc-FPR1) + TPR2*(FPR2-cc) ) / (FPR2-FPR1)
		FPRcut=c(FPR[1:FPRlcc],cc)
		TPRcut=c(TPR[1:FPRlcc],TPRcc)
	}
	return(sum( diff(FPRcut) * ( TPRcut[-1] + head(TPRcut,-1) ) )/2)
}

AUC <- function(prob, outcome, cutoff = 1, ROC.plot = FALSE)
{
	NN <- length(outcome)
	if(length(prob)!=NN) stop("prob and binary should be the same length")
	if((max(prob)>1) | (min(prob)<0)) stop("prob values should be in [0,1]")

	# remove missing values, which cannot contribute to AUC calculation
	idx.retain <- (!is.na(outcome)) & (!is.na(prob))
	prob <- prob[idx.retain]
	outcome  <-  outcome[idx.retain]

	# sort binary outcome by order of predicted probability
	ord <- order(prob, decreasing=T)
	prob <- prob[ord]
	outcome <- outcome[ord]

	CP <- sum(outcome)    # condition positive
	CN <- NN-sum(outcome) # condition negative
	TP <- cumsum(outcome) # true positive
	FP <- (1:NN) - cumsum(outcome) # false Positve
	TPR <- TP/CP		 # True positive rate (sensitivity)
	FPR <- FP/CN		 # False positive rate (1- specificity)
	if(ROC.plot) {plot(FPR, TPR); abline(v=cutoff)}
	auc <- trapezoid(FPR, TPR, cutoff)

	return(auc)
}
