

scores <- function(predictions, testlabels) {
  # F_1 = 2*(TP) / {2*(TP) + FN + FP}
  pos.indices <- which(testlabels %in% 1)
  neg.indices <- which(testlabels %in% 0)
  
  TP = sum(predictions[pos.indices] == 1)  ### out of n
  TP_rate = 
  TN = sum(predictions[neg.indices] == 0)
  FN = sum(predictions[pos.indices] == 0)
  FP = sum(predictions[neg.indices] == 1)
  FP_rate = 
  precision <- TP/(TP + FP)
  recall <- TP/(TP + FN)
  F1 <- 2 * (precision*recall/(precision + recall))  
  lrplus <- TP/FP
  lrminus <- FN/TN
  dor <- lrplus/lrminus
  
  return(c(TP_rate, FP_rate, F1, lrplus, lrminus, dor))
}


#####  1.)  Covariance Function     ########################################
##
##
##
##
covariance <- function(level, Predictors, noise)
{
  
  
  if (level=='low') {cov.level<-.2
                     
  } else if (level == 'lower') { cov.level<- 0.1
                                 
  } else if (level == 'zero') { cov.level<- 0
                                
  } else if (level == 'low/med') { cov.level<- .3                              
                                   
  } else if (level == 'med') { cov.level<-.5
                               
  } else if (level == 'high'){ cov.level<-.8
                               
  } else if (level == 'higher'){ cov.level<-.9                             
                                 
  }
  
  
  
  covmat<-matrix(1, Predictors, Predictors)
  
  covmat[upper.tri(covmat)]<-cov.level
  
  covmat<-covmat+ matrix(rnorm(Predictors^2, 0 , noise), Predictors, Predictors)
  
  covmat[!upper.tri(covmat)]<-0
  
  covmat<-covmat+t(covmat)
  
  diag(covmat)<-1
  
  
  return(covmat)
  
}




###           2.) Sampling Function       #################################
##
##   Sample artificial data for given environment 
##  function 'sampling' samples binary x and binary y
##
##
##
sampling <- function(Predictors, N, CovMat, noise, betas)
{
  
  ##   Sampling is a function where:
  ##  	for each trial, 
  ##    	a. sample values of x based on the covariance matrix of step 1
  ##    	b. sample y based on the regression weights from step 2
  ##
  ##	Arguments:
  ##
  ##	betas			- sampled regression weights (Gaussian/Laplace)
  ##	Predictors 		- no. of predictors
  ##	N	    	 	- sample size 
  ##	CovMat		- covariance matrix
  ##	seed			- current seed (for current simulation)
  ##	Noise 		- noise level (small, medium, large)
  ##
  ##	Return values: the matrix of binary x's (NxPredictors) and classified Ys
  ##	vector (Nx1 vector) as a combined data frame
  
  ##	mvrnorm simulates N samples from a multivariate normal distribution 
  ##	with the covariances between Predictors according to Sigma matrix
  ##	mu is a vector of means=0 for the number of Predictors
  ##	CovMat uses covariance matrix with 1's in diagonal, and covariances between 0 and 1
  
  ## 	setting the seed for mvrnorm() gives the same binary data marix
  ##	'bdata' each time,but be careful not to use the same as for betas!
  
  #   seed_int <- round(runif(1, min=1, max = 200)) # random seed
  #   set.seed(seed_int)
  
  mu <- rep(0,Predictors) # all means are zero
  
  
  ##    Create a differences matrix beween cue matrix 1 and 2 (option A and B)
  ##    pairs of items are collapsed down into one matrix
  
  # 1. sample cue matrix 1 (option A) with covaraiance sigma = CovMat
  
  mdata1 <- mvrnorm(n=N, mu = mu, Sigma = CovMat, empirical = TRUE) 
  
  cat("Covariance matrix used for predictors:\n\n")  
  print(var(mdata1)) # print the covariances of the sampled data
  
  ## generate binary data from multivariate data 
  
  bdata1 <-(mdata1 > 0)*1
  
  # 2. sample cue matrix 2 (option B) with covariance sigma = CovMat
  
  mdata2 <- mvrnorm(n=N, mu = mu, Sigma = CovMat, empirical = TRUE) 
  
  ## generate binary data from multivariate data 
  
  bdata2 <-(mdata2 > 0)*1
  
  
  # 3. differences matrix [matrix 1 - matrix 2]
  
  bdata_diff <- bdata1 - bdata2
  
  
  
  ## 	4. matrix multiplication
  
  
  Xw <- bdata_diff %*% betas
  
  
  ##	2. add Gaussian noise  
  
  
  Ymulti <- Xw + noise
  
  
  ###TALLYING ENVIRONMENT #################################
  ## thresholding Ymulti at 0 for differences matrix bdata_diff (Gaussian)
  
  Ymulti[Ymulti > 0]  <- 1
  Ymulti[Ymulti < 0]  <- - 1
  
  
  # when tallying would guess
  for (j in 1:length(Ymulti)){  
    if (Ymulti[j]==0)
    {
      Ymulti[j] <- sign(rnorm(1))
    }
  }   ## all +1 and -1 now!
  
  
  y <- Ymulti
  #y <- (Ymulti > 0) * 1
  
  
  
  environment <- data.frame(bdata_diff, dependent=y)
  
  return(environment)
  
}





########## 	        3.) Function cv.indexing                        ################################
##
##	cross-validation partitioning for current environment dataset,
##	takes seed as input to be able to reproduce results,
##	However, can use different seeds here to get different cross 
##	validation results for the same data set 

## Generates 'k' sets of indices, where each set represents a split
## of a data set of size 'N' into two disjoint data sets (sometimes called a
## 'folder').
##  
## Within each split, part of the data is called a 'test' data set, and the
## remaining data points form the 'training' data. None of the 'test' sets,
## overlap, and their union corresponds to the complete data set.
##
## This function returns a k-by-n binary matrix 'cv'. Each row in 'cv' 
## represents a split. If cv[i, j] == 1, then data point 'j' is in the 
## training set of partition 'i'. If cv[i, j] == 0, then data point 'j' is 
## in the test set of partition 'i'.
##
## The number of test points in the last partition will be larger than in 
## the other partitions if 'N' is not a multiple of 'k'.

cv.indexing <- function(k, N, percent)
{
  ## 	Arguments:
  ##	k 	number of partitions (set above)
  ##	N	size of the current data set (set above)
  ##	seed 	seed used for random partitions (set above) 
  
  #   if (k >= N)
  #     stop('Number of partitions should be less than number of data points!')   
  #   
  
  #set.seed(seed) # For reproducibility of results, only
  # dont use it now for the perm function!
  
  test.set<- round(percent * N)
  # the size of a test set is 50% of N --> abc data sets!
  
  
  cv <- matrix(rep(1, k * N), nrow = k) # Allocate output matrix
  
  
  for (i in 1:k){
    
    # generate new permuation for each row i 
    perm <- sample(1:N, N)
    
    # define the test set items from columns in perm, to be marked by 0 
    # while the training set items stay 1
    cv[i,perm[1:test.set]] <- 0
    
    
  }
  
  return(cv)                              
}

####################    3.)  Logistic Regression    ########################
###     as traditional linear regression comparison model
##
##
##

logistic.regression <- function(paired_data, cv, cv2, y.pos) 
{
  
  k <- nrow(cv) 
  pred.logistic <- vector('list', k)
  
  for (i in 1:k){
    
    trainset <- paired_data[cv2[i,]==1, ]
    variables <- names(trainset) # of training set only!
    var <- "factor(dependent, levels= c(0,1))"

     #----------------     Fitting the training data  ----------------------------#########################

    logistic.fit <- glm(paste(var," ~ ", paste(variables[-y.pos], collapse= "+"), paste(" -1")), 
                        family=binomial, data=trainset)
    summary(logistic.fit)
    
    #-----------------  Predicting the Test set ------------------------------------------##########
    
    testset <- paired_data[cv2[i,]==0, ]
    
    # type "class" does not work for predict.glm, only response or link (log odds)
    # type = response gives predicted probabilities of y=1 (compared to y=-1)
    log.predict <- as.vector(predict(logistic.fit,newdata=testset[ ,-y.pos], type = "response"))
    
    # threshold the probability > 0.5 to y=1 otherwise y=-1, has to be done this way 
    # otherwise some values will stay 0.5 and not classified
    ## is equivalent to ridge and lasso threshold
    log.predict <- (log.predict > 0.500000000)*1
    # this time only 0 and 1 necessary
    
    
    pred.logistic[[i]] <- log.predict
    
    
  }
  
  return(pred.logistic)   
  
}


####################   4.)  Logistic Regression with threshold parameter   ########################
###     as traditional linear regression comparison model
##
##
##

# logistic.regression <- function(paired_data, cv, cv2, y.pos, Predictors) 
# {
# 
#   #-------------- Find the best threshold cut-off for logistic regression  -------------------
#   k <- nrow(cv)  
#   # threshold parameter values to test
#   threshold <- seq(0.1,1,0.05)   
#   acc.t <- vector('numeric', length(threshold))
#   # t runs through all possible threshold values    
#   for (t in 1:length(threshold)){   
#     
#     thres.acc <- vector('numeric', k)
#     
#     for (i in 1:k){   
#       
#       trainset <- paired_data[cv[i,]==1, ]
#       testset <-paired_data[cv[i,]==0, ]
#       actual_labels <- testset[ ,ncol(paired_data)] # optimizing performance in the test set
#         
#       ### Fitting the training set ################################################
#       variables <- names(trainset) # of training set only!
#       fmla <- as.formula(paste(variables[y.pos]," ~ ", paste(variables[-y.pos], collapse= "+"), paste(" -1")))     
#       logistic.fit <- glm(fmla, family=binomial, data=trainset)    
#       log.predict <- as.vector(predict(logistic.fit,newdata=testset, type = "response"))
#       
#       # threshold the probability at different thresholds t (from 0.1 to 1) and if the probability > threshold, 
#       # y=1 otherwise y=0 (1= present, 0=absent)    
#       log.predict <- (log.predict > threshold[t])*1
#       # get the threshold accuracy for this testset
#       thres.acc[i]   <- regression.graph(log.predict, actual_labels)  
#      
#     }
#      # get the overall accuracy for each threshold (across trainsets)
#     acc.t[t] <- sum(thres.acc)/k
#     
#   }
#   
#   ### Optimal threshold is:   
#   max.pos  <- which.is.max(acc.t)
#   optimal.t <- threshold[max.pos]
#   
#   
#   ###########-------------- Use the optimal threshold to predict 20 new cross-validation (cv2) sets -----------
#   m <- nrow(cv2)     
#   pred.logistic <- vector('list', m)
#   
#   for (j in 1:m){
#       
#     trainset2 <-  paired_data[cv2[j, ]==1, ]
#     testset2 <-  paired_data[cv2[j,]==0, ] 
#      
#     # ------------- Fitting the training set ---------------------------------------------------------###########
#     variables <- names(trainset2) # of training set only!
#     fmla <- as.formula(paste(variables[y.pos]," ~ ", paste(variables[-y.pos], collapse= "+"), paste(" -1")))   
#     
#     
#     actual.fit <- glm(fmla, family=binomial, data=trainset2)  
#     summary(actual.fit)
#     
#     #---------- Deriving the predictions for the new test sets and threshold them at optimal threshold -----------------------
#     final.predict <- as.vector(predict(actual.fit,newdata=testset2, type = "response"))    
#     # threshold the probability at optimal thresholds optimal.t
#     # if the probability > threshold,  y=1 otherwise y=0 (1= present, 0=absent) 
#     final.predict <- (final.predict > optimal.t)*1
#   
#     # entre each final.predict vector into list element i
#     pred.logistic[[j]] <-  final.predict
#     
#   } # j loop
#   
#   
#   AllList <- list(predictions=pred.logistic, optimal_threshold = optimal.t)
#     
#   return(AllList)  
#   
# }
#   


####################    4.)  support Vector Machine (Gaussian Kernel)   ########################
##
##  to test for interactions
##

# 
# gamma <- 2^(seq(-15,3,2))
# C <- 2^(seq(-5,15,2))
# y.pos <- ncol(paired_data)

svm.method <- function(paired_data, cv, cv2, y.pos,gamma, Cost) 
{
  

  k <- nrow(cv)     
  all.accuracies <- vector('list', k)# list that can contain vectors with performances of all parameter 
 
  ## Create all possible combinations of gamma and C: Grid Search
  # a data frame containing one row for each possible parmater combination
  parameters <- expand.grid(Cost=Cost, gamma=gamma)

  
  for (i in 1:k){
    
    ## Location in the loop to the prompt
    cat(paste("\n\nTEST SET NUMBER  = ", i, "\n\n"))
    trainset <- paired_data[cv[i,]==1, ]
    testset <- paired_data[cv[i,]==0, ]
    
    svm.pred <- matrix(nrow=nrow(parameters),ncol=nrow(testset))  # 110 combinations x nrow(testset)   
 
    for (p in 1:nrow(parameters))
    {
      
      #-------------------  FITTING SVM with one parameter combination (p) at a time --------------------------- ######################################    
      param_tune <- parameters[p, ]
 
      # scale = TRUE      
#       ## this is not classification but eps -regression!!!
#       svm.model <- svm(factor(dependent,levels= c(0,1)) ~ ., data = trainset, 
#                         cost = param_tune$Cost, gamma =param_tune$gamma)


      # dependent is now already a factor (for C-classification!) 
      svm.model <- svm(dependent ~ ., data = trainset, kernel = "radial", 
                       cost = param_tune$Cost, gamma =param_tune$gamma)
      
#       ------------------------- PREDICTING TESTSET WITH SVM MODEL (parameters p) -------------------------
#       Warning: To convert factors to numeric or integer, first convert to character. 
#       Converting factors directly to numeric or integer data can lead to unwanted outcomes. 

      svm.char <- as.character(predict(svm.model, testset[ ,-y.pos]))
      svm.pred[p, ] <- as.numeric(svm.char) # characters have to be converted via "as.numeric" instead of "as.integer"
           
    } # parameters loop p
    
    #####----------------- Get accuracy for EACH PARAMETER COMBINATION (for the CURRENT TEST SET i) ----------------####      
    # compare each vector of predictions of all parameters p to actual lables (factor is ok here!)
    
    test.labels <- testset[ ,ncol(paired_data)]
#svm.predict is referring to the p x n(testset) matrix
    all.accuracies[[i]] <- pred.graph(svm.pred, test.labels)
    # contains each one list i with vector for each testset i $raw.acc
     
  } # k loop

  
  ## overall accuracy across k test sets, per penalty
  accuracies <- vector("numeric",k)
  acc <- vector('numeric', nrow(parameters))
  
  for (l in 1:nrow(parameters)){   
    for (i in 1:k) {
      accuracies[i]   <- all.accuracies[[i]]$raw.acc[l]     
    }    
    acc[l]  <- sum(accuracies)/k        
    ## acc contains accurcay of parameter combination p across all train-test set separations!
  }
  
  
  ####----- OPTIMAL Parameter Combination across all these 20 test sets-----------------------------------------
  
  # get the maximum performance penalty
  max_parameters <- which.is.max(acc) # location of optimal penalty
  optimal_parameters <- parameters[max_parameters, ] 
  
  # get the accuracy(performance) at that max penalty
  max_acc <- acc[max_parameters]  
  #----------------------------------------------------------------------------------------------
  
  
  ###########-------------- Use the optimal parameters to predict 20 new cross-validation sets -----------
  
  m <- nrow(cv2)     
  pred.svm <- vector('list', m)
  prediction.accuracies <- vector('numeric', m)# list that can contain vectors with performances of all parameter 

  for (j in 1:m){
    
    trainset2 <- paired_data[cv2[j,]==1, ]
    testset2 <- paired_data[cv2[j,]==0, ]
    
    # fit with optimal parameter combination from above and new training data
    svm.fit <- svm(dependent ~ ., data = trainset2, kernel = "radial", 
                     cost = optimal_parameters$Cost, gamma = optimal_parameters$gamma)
 
    #       ------------------------- PREDICTING TESTSET WITH SVM MODEL (parameters p) -------------------------
    #       Warning: To convert factors to numeric or integer, first convert to character. 
    #       Converting factors directly to numeric or integer data can lead to unwanted outcomes.  
    svm.charpred <- as.character(predict(svm.fit, testset2[ ,-y.pos]))
    pred.svm[[j]] <- as.numeric(svm.charpred)                  
    
#      # to get a peak of the accuracies and compare with the one in the main script!
#       test.labels <- testset2[ ,ncol(paired_data)]
#       prediction.accuracies[j] <- regression.graph(pred.svm[[j]], test.labels) 

  } # cv2 loop (j)
  
  # test
  #overall <- sum(prediction.accuracies)/k


  ##--------------- to get the support vector coefficients fit the whole dataset once 
  # ?????? or is there anyway to average them across cv sets?
  # fit paired_data with optimal parameters gamma and COst 
  whole.fit <- svm(dependent ~ ., data = paired_data, 
                   cost = optimal_parameters$Cost, gamma = optimal_parameters$gamma) 
  sv_coeff <- whole.fit$coefs
  ## -----------------------------------------------------------------------------
  

  AllList <- list(svm_predictions=pred.svm, optimal_parameters = optimal_parameters,
                  sv_coeff = sv_coeff)
  
  return(AllList)

}







####################    3.)  Logistic Regression Manual (OLS by hand)   ########################
##
##
##

logistic.manual <- function(paired_data, cv, y.pos, Predictors) 
{
  
  k <- nrow(cv) 
  pred.logistic <- vector('list', k)
  
  for (i in 1:k){
    
    trainset <- paired_data[cv[i,]==1, ]
    
    #### Estimate the LINEAR regression weights (betas) with OLS (least squares)
    ##### after binarizing; stricly y is binary -1/+1 so only Logistic Regression OLS estimate is legite 
    X <- as.matrix(trainset[  ,1:Predictors])
    Y <- as.matrix(trainset[ ,ncol(trainset)])
    
    print(is.non.singular.matrix(t(X)%*% X))
    
    if (is.non.singular.matrix(t(X)%*% X) == 'FALSE'){
      
      weights <- rep(0, Predictors)
      
    } else{
      
      weights <- as.vector(solve(t(X)%*%X) %*% (t(X)%*%Y))  # OLS estimation of parameters
      # roughly noncompensatory weights should be 0.8,0.4,0.2            
    }
    
    #-----------------  Predicting the Test set ------------------------------------------##########
    testset <- paired_data[cv[i,]==0, ]
    test <- as.matrix(testset[ ,1:Predictors])    
    logpredict <- as.vector(t(weights %*% t(test))) 
    
    # thresholding the matrix multiplications to y predictions 1/-1
    logpredict[logpredict > 0] <- 1
    logpredict[logpredict < 0] <- -1
    
    # put vector of predictions into list element i
    pred.logistic[[i]] <- logpredict
    
  }
  
  return(pred.logistic)   
}





#############       4.) Regression Prediction Graph (accuracies)########################

regression.graph <- function(predictions.regr, test.labels)
{
  
  #predictions.regr <- predictions.scale[[i]]
  #predictions.regr <- predictions.regression[[1]]
  
  # predictions.regression      is a vector from a list containing the 0/1 predictions
  # test.labels                 is a vector of the correct outcomes from the current testset
  
  n <- length(predictions.regr)
  m <- length(test.labels)
  all.equal(n,m)
  
  prop_regr <- sum(test.labels==predictions.regr)/m
  
  return(prop_regr)
  
}



########    5.) Regular Ridge Function (Gaussian, non truncated, non Bayesian) ###############################################
#        
# for testing only
#y.pos <- ncol(paired_data)

regular.ridge <- function(paired_data, cv, cv2, y.pos, penalty) 
{ 
  
  k <- nrow(cv)     
  all.accuracies <- vector('list', k)
  
  for (i in 1:k){
    
    ## Location in the loop to the prompt
    cat(paste("\n\nTEST SET NUMBER  = ", i, "\n\n"))
    
    trainset <- paired_data[cv[i,]==1, ]
    testset <- paired_data[cv[i,]==0, ]
    
    # Design matrix X of TEST SET
    X <- as.matrix(trainset[-y.pos]) # has to be a matrix   
    # dependent y is a factor(!) with levels 0/1
    y <- factor(trainset[ ,y.pos], levels= c(0,1))
    
    log.predict <- matrix(nrow=length(penalty),ncol=nrow(testset))     # 33xnrow(testset)     
    
    for (p in 1:length(penalty))
    {
           
      #-------------------  FITTING RIDGE REGRESSION  with 1 penalty at a time --------------------------- ######################################    
      # alpha = 0 indicates ridge
      # family = binomial indicates logistic regression
      # intercept = FALSE to stay consistent with the glm fit for logistic regression
      
      logRidge.fit <- glmnet(X, y, family="binomial", alpha = 0, lambda=penalty[p],intercept=FALSE)
      #logRidge.fit$beta


      #-----------------  PREDICTING with that same penalty (predicted probabilities) --------------------------################          
      # get test data 
      testmatrix <- as.matrix(testset[-y.pos])
      # naming the penalty is not necessary again, it is stored in logRidge.fit
      # type= "class" corresponds to "response" probabilities         
      #log.predict[p, ] <- as.vector(predict(logRidge.fit, newx = testmatrix,type="class"), mode= "numeric")        
      
      temp <- as.vector(predict(logRidge.fit, newx = testmatrix,type="response"), 
                                    mode= "numeric")  
      
      # threshold equvialently to linear regression:
      #the probability > 0.5 to y=1 otherwise y=-1, has to be done this way 
      temp <- (temp > 0.500000000)*1
      # this time temp only contains 0 and 1 (physical aggression present/absent)
      
      
      log.predict[p, ] <- temp
      
      
    } # penalty loop p
       
    #####--------------------- Get accuracy score of EACH PENALTY (for the CURRENT TEST SET i) ----------------####      
    # returns a vector of accuracies for all 33 penalties
    # lop.predict is matrix with penalty x nrow(testset) dimensions
    
    test.labels <- testset[ ,ncol(paired_data)]
    all.accuracies[[i]] <- pred.graph(log.predict, test.labels)
    
  } # k loop


  ## overall accuracy across k test sets, per penalty
  accuracies <- vector("numeric",k)
  acc <- vector('numeric', length(penalty))
  
  for (l in 1:length(penalty)){   
    for (i in 1:k) {
      accuracies[i]   <- all.accuracies[[i]]$raw.acc[l]     
    }    
    acc[l]  <- sum(accuracies)/k        
    ## acc contains averaged accuracies per penalty!
  }
  
  
  ####----- OPTIMAL PENALTY across all these 20 test sets-----------------------------------------
  
  # get the maximum performance penalty
  max_p <- which.is.max(acc) # location of optimal penalty
  optimal_p <- penalty[max_p] 
  
  # get the accuracy(performance) at that max penalty
  max_acc <- acc[max_p]  
#----------------------------------------------------------------------------------------------


###########-------------- Use the optimal penalty to predict 20 new cross-validation sets -----------

  m <- nrow(cv2)     
  pred.logridge <- vector('list', m)

  for (j in 1:m){

    trainset2 <- paired_data[cv2[j,]==1, ]
    testset2 <- paired_data[cv2[j,]==0, ]

    # Design matrix X of TEST SET
    X <- as.matrix(trainset2[-y.pos]) # has to be a matrix   
    y <- factor(trainset2[ ,y.pos], levels= c(0,1))

    # fit with optimal penalty and new training data
    actual.fit <- glmnet(X, y, family="binomial", alpha = 0, lambda=optimal_p,intercept=FALSE)    
    
#-----------------  PREDICTING with that same penalty (predicted probabilities) --------------------------################          
    
    test2 <- as.matrix(testset2[-y.pos]) # test data from new testset2    
    pred.logridge[[j]] <- as.vector(predict(actual.fit, newx = test2,type="class"), 
                                  mode= "numeric")  
   
  } # end of cv2 loop (j)



##--------------- to get the optimal regularized coefficients fit the whole dataset once 

X2 <- as.matrix(paired_data[-y.pos])
y2 <- factor(paired_data[ ,y.pos], levels= c(0,1))
# fit with optimal penalty and new training data
whole.fit <- glmnet(X2, y2, family="binomial", alpha = 0, lambda=optimal_p,intercept=FALSE)    
optimal_coeff <- whole.fit$beta

## -----------------------------------------------------------------------------

  
  AllList <- list(optim_predictions=pred.logridge, optimal_p = optimal_p,
                  optimal_coeff = optimal_coeff)
  
  
  return(AllList)
  
}  




#####     6.) Pred.Graph Function (Ridge accuracies)    #########################
##
##
##
pred.graph <- function(predictions, test.labels){
  
  # predictions <- log.predict        ### is matrix with penalty x nrow(testset) dimensions
  # predictions <- svm.pred
  
  #   get the number of test predictions 
  n <- ncol(predictions)
  
  raw.acc <- vector("numeric",nrow(predictions))
  
  for (l in 1:nrow(predictions)){
    
    raw.acc[l] <- sum(predictions[l, ] == test.labels)/n
    
  }
  
  # put the vector raw_acc into a list each time, and return it	
  
  propList <- list(raw.acc=raw.acc)
  return(propList)
  
}



############ 7.) LASSO PREDICTION FUNCTION ####################################
##
##
##

regular.lasso <- function(paired_data, cv, cv2, y.pos, penalty) 
{
  
  k <- nrow(cv) 
  all.accuracies <- vector('list', k)
  
  for (i in 1:k){
    
    ## Location in the loop to the prompt
    cat(paste("\n\nTEST SET NUMBER  = ", i, "\n\n"))
    
    trainset <- paired_data[cv[i,]==1, ]
    testset <- paired_data[cv[i,]==0, ]
    
    # Design matrix X of TEST SET
    X <- as.matrix(trainset[-y.pos]) # has to be a matrix   
    # dependent y is a factor(!) with levels -1/1
    y <- factor(trainset[ ,y.pos], levels= c(0,1))
    
    log.predict <- matrix(nrow=length(penalty),ncol=nrow(testset))     # 33xnrow(testset)     
   
    for (p in 1:length(penalty))
    {
    
     #-------------------  FITTING lasso REGRESSION  with 1 penalty at a time --------------------------- ######################################    
      # alpha = 1 indicates lassom, family = binomial indicates logistic regression
      
      loglasso.fit <- glmnet(X, y, family="binomial", alpha = 1, lambda=penalty[p],intercept=FALSE)
      
      #-----------------  PREDICTING with that same penalty (predicted probabilities) --------------------------################          
      # get test data 
      testmatrix <- as.matrix(testset[-y.pos])
      # naming the penalty is not necessary again, it is stored in loglasso.fit
      # type= "class" corresponds to "response" probabilities         
      temp <- as.vector(predict(loglasso.fit, newx = testmatrix,type="response"), 
                                    mode= "numeric")        
     
     # threshold equvialently to linear regression:
     #the probability > 0.5 to y=1 otherwise y=-1, has to be done this way 
     temp <- (temp > 0.500000000)*1
     # this time temp only contains 0 and 1 (physical aggression present/absent)
     
     
     log.predict[p, ] <- temp
        
      
    } # penalty loop p
    

    #####--------------------- Get accuracy score of EACH PENALTY (for the CURRENT TEST SET i) ----------------####      
    # returns a vector of accuracies for all 33 penalties
    
    test.labels <- testset[ ,ncol(paired_data)]
    all.accuracies[[i]] <- pred.graph(log.predict, test.labels)
    
  } # k loop
  
    
    ## overall accuracy across k test sets, per penalty
    accuracies <- vector("numeric",k)
    acc.lasso <- vector('numeric', length(penalty))
    
    for (l in 1:length(penalty)){
      for (i in 1:k) {
        accuracies[i]   <- all.accuracies[[i]]$raw.acc[l]     
      }    
      acc.lasso[l]  <- sum(accuracies)/k        
      ## acc contains averaged accuracies per penalty!
    }
  
  
  ####----- OPTIMAL PENALTY -----------------------------------------
  
  # get the maximum performance penalty
  max_p <- which.is.max(acc.lasso)
  optimal_p <- penalty[max_p] 


    ###########-------------- Use the optimal penalty to predict 20 new cross-validation sets -----------
    
    m <- nrow(cv2)     
    pred.loglasso <- vector('list', k) 
    
    for (j in 1:m){
      
      trainset2 <- paired_data[cv2[j,]==1, ]
      testset2 <- paired_data[cv2[j,]==0, ]
      
      # Design matrix X of TEST SET
      X <- as.matrix(trainset2[-y.pos]) # has to be a matrix   
      y <- factor(trainset2[ ,y.pos], levels= c(0,1))
      
      # fit with optimal penalty and new training data
      actual.fit <- glmnet(X, y, family="binomial", alpha = 1, lambda=optimal_p,intercept=FALSE)      
      
      #-----------------  PREDICTING with that same penalty (predicted probabilities) --------------------------################          
      
      test2 <- as.matrix(testset2[-y.pos]) # test data from new testset2    
      pred.loglasso[[j]] <- as.vector(predict(actual.fit, newx = test2,type="class"), 
                                      mode= "numeric")  
      
    } # end of cv2 loop (j)
    
    
    ##--------------- to get the optimal lasso regularized coefficients fit the whole dataset once 
    X2 <- as.matrix(paired_data[-y.pos])
    y2 <- factor(paired_data[ ,y.pos], levels= c(0,1))
    # fit with optimal penalty and new training data
    whole.fit <- glmnet(X2, y2, family="binomial", alpha = 1, lambda=optimal_p,intercept=FALSE)    
    optimal_coeff <- whole.fit$beta   
    ## -----------------------------------------------------------------------------
    
    
    AllList <- list(optim_predictions=pred.loglasso, optimal_p = optimal_p,
                    optimal_coeff = optimal_coeff)
    
  
  return(AllList)
  
}


#####     8.) Lasso.Graph Function (lasso accuracies)    #########################
##
##
##
# lasso.graph <- function(predictions, test.labels){
#   
#   
#   # predictions <- predictions.only[[i]] is a character matrix still!!  
#   
#   predictions <- as.vector(predictions, "numeric") # convert to numeric vector
#   
#   n <- length(predictions) # vector
#   m <- length(test.labels)
#   all.equal(m,n)
#   
#   # is a scalar of %
#   lasso.acc <- sum(predictions == test.labels)/n
#   
#   return(lasso.acc)
#   
#   
# }


############ 7.) SIMPLE SCALE (THRESHOLD MODEL) ####################################
##
##
##   optimizing the threshold paramter to perform optimal in the test set or training set
##   Has to be set below!
 
# y.pos <- ncol(scale_data)

simple.scale <- function(scale_data, cv, cv2, y.pos, Predictors) 
{
 
  
  k <- nrow(cv)  

  # get the minimum and maximim value in the sample to get range of threshold to test for
  min <- min(scale_data$hcr_total) # minimum value in the sample
  max <- max(scale_data$hcr_total) # maximum value in the sample
  
  threshold <- 1: (2*Predictors) # maximum score possible with amount of predictors! 
  
  
  acc.t <- vector('numeric', length(threshold))
  
  # t runs through all possible threshold values    
  for (t in 1:length(threshold)){   
    
    thres.acc <- vector('numeric', k)
    
    for (i in 1:k){   
      
    trainset <- scale_data[cv[i,]==1, ]
    testset <- scale_data[cv[i,]==0, ] 
    #score <- trainset[ ,1] # optimizing perfromance in the train set
    score <- testset[ ,1] # optimizing performance in the test set
    #actual_labels <- trainset[ ,2] # optimizing performance in the train set
    actual_labels <- testset[ ,2] # optimizing performance in the test set

      score_predict   <- (score > threshold[t]) * 1          
      thres.acc[i]   <- regression.graph(score_predict, actual_labels)  
   
    }
    
    
    acc.t[t] <- sum(thres.acc)/k
       
  }
  
    max.pos  <- which.is.max(acc.t)
    optimal.t <- threshold[max.pos]
    
    
   
  ###########-------------- Use the optimal threshold to predict 20 new cross-validation sets -----------
      
      m <- nrow(cv2)     
      pred.scale <- vector('list', m)
      
      for (j in 1:m){
           
        
        trainset2 <- scale_data[cv2[j]==1, ]
        testset2 <- scale_data[cv2[j,]==0, ] 
        #score <- trainset[ ,1] # optimizing perfromance in the train set
        score_test <- testset2[ ,1] # optimizing performance in the test set
  
        ## get predictions: Thresholding the test set score at the opitmal threshold optimal.t
        test.predictions <- (score_test > optimal.t)*1
        
        # entre each test.predictions vector into list element i
        pred.scale[[j]] <-  test.predictions
    

    } # j loop
  
  ## RETURN THE OPTIMAL THRESHOLD!
  return.stuff <- list(pred.scale = pred.scale, optimal.threshold = optimal.t)
  
  return(return.stuff)  
  
}
  


#####     9.) Accordance Function (Heuristic and Regularized)    #########################
##
##
##
accordance <- function(predictions.heuristic, predictions.regularized)
{
  
  # predictions.regularized <- ridge_matrix[i, ]   is already a numeric vector row by row (i)
  # predictions.regularized <- lasso_matrix[i, ]   is already a numeric vector row by row (i)
  
  
  
  n <- length(predictions.regularized) # vector
  m <- length(predictions.heuristic)
  all.equal(m,n)
  
  # is a scalar of %
  acc_rate <- sum(predictions.regularized == predictions.heuristic)/n
  
  
  
  return(acc_rate)
  
}




########## 9.) Tallying Learning Function ##############################


#y.pos <- ncol(paired_data)

tallying.learning <- function(paired_data, cv2, y.pos, Predictors) 
{
  
  #	number k of different test sets
  k <- nrow(cv2) 
  
  ## 	initiate empty lists
  pred.tallying <- vector('list', k)
  
  for (i in 1:k){
    
    ####  Fitting Training Data (get unit weights) ###################################
    trainset <- paired_data[cv2[i,]==1, ]
    # use the ecological cue validities, but make them +1 if the validity is above chance (v > 0.5)
    # and -1 if it is below chance
    
    cue_validities <- vector('numeric', Predictors)
    for (c in 1:Predictors){
      # estimate the ecological cue validity of each cue as v = R/(R+W)
      if (sum(trainset[,c]==trainset[ ,ncol(trainset)]) == 0) { cue_validities[c] <- 0 
      } else  
        cue_validities[c] <- sum(trainset[,c]==trainset[ ,ncol(trainset)])/(sum(trainset[,c]==1)+sum(trainset[,c]==-1))              
    }
    
    
    # these unit weights would be different if regression weights would be used instead of cue validities (1, 0.5. 0.25)
    
    
    # derive +1/-1 unitweights: make them -1 if the weight was below 0.5
    unitweights <- (cue_validities > 0.5)*1
    unitweights[(unitweights == 0)] <- -1    
    
    ############ Prediction Tallying (Testdata) ###########################################    
    testset <- paired_data[cv2[i,]==0, ]
    test <- as.matrix(testset[ ,-y.pos]) #  should be 1702x9 double matrix, without dependent
    
    # matrix multiplication of unit weights and testdata, if > 0, then +1 (A), if < 0, then B(-1)
    tl_predict <- sign(test %*% unitweights)
    
    ## if the sum of all cue differences is still 0, then guess, different each time!!
    for (j in 1:length(tl_predict)){  
      if ( tl_predict[j]==0)
      {
        tl_predict[j] <- sign(rnorm(1))
      }
    }   ## all +1 and -1 now! CHNAGE TO 0 and 1
    
    
    # feed vector for one testset into predictions list
    pred.tallying[[i]] <-  tl_predict 
    
  }
  
  return(pred.tallying)
}




#####         9.) Tallying.graph (accuracies) Function    ######################################
##

tallying.graph <- function(predictions.tallying, test.labels)
{
  
  ## function provides a summary of the performance (accuracy) of the tallying 
  ## heuristic
  n <- length(predictions.tallying)
  m <- length(test.labels)
  all.equal(n,m)
  
  prop <- sum(test.labels==predictions.tallying)/m
  
  return(prop) # a single number 
  
}



#### 10. TTB Learning Function ##########################################################################

# data <- paired_data
#  y.pos <- ncol(paired_data)

ttb.predictions <- function(data, cv, y.pos, Predictors)
{
  ##   This function is called to derive the predictions for Take-the-Best 
  ## Arguments:
  ##
  ## dataset 		dataset is data frame of current environment with n rows = data points
  ## cv 	  	  a binary matrix with k x n; k rows correspond to 
  ##			      different test and training set separations, 
  ##			      n columns correspond to the n data points in the dataset
  ## y.pos 		  integer indicating which of the columns in dataset 
  ##			      contains the binary dependant of the data points (0/1)
  
  #  number k of different test sets
  k <- nrow(cv) 
  
  ## 	initiate empty lists
  pred.ttb <- vector('list', k) # because we have several test sets here
  
  #-------------------------- get the cue validities and order (Trainset) ---------------------------------------
  for (i in 1:k){
    
    trainset <- data[cv[i,]==1, ]
    
    cue_validities <- vector('numeric', Predictors)
    for (c in 1:Predictors){
      # estimate the ecological cue validity of each cue as v = R/(R+W)
      if (sum(trainset[,c]==trainset[ ,ncol(trainset)]) == 0) { cue_validities[c] <- 0 
      } else  
        cue_validities[c] <- sum(trainset[,c]==trainset[ ,ncol(trainset)])/(sum(trainset[,c]==1)+sum(trainset[,c]==-1))             
    }
    
    
    # cue validities are all > 0, between 0-1
    cue_order <- order(cue_validities, decreasing = TRUE)
    
    ############ PREDICTING TESTDATA WITH TTB (PREDICTION) ###########################################      
    testset <- data[cv[i,]==0, -y.pos]
    ttb_predict <- vector("numeric",nrow(testset))   
    
    # TTB Search Algorithm
    for(r in 1:nrow(testset)){
      v <- 1
      # circulates through elements 1:5 of cue_order in order of abs. vailidity 
      while (testset[r,cue_order[v]]== 0){            
        
        if (v == Predictors){          
          ttb_predict[r] <-  sign(rnorm(1)) # if no cue discriminates, guess between A(+1) and B(-1)
          break                
        }       
        v <- v + 1  # go to the next element in the cue_order  
      } # while loop only breaks when testset ~= 0 
      
      ## check if the cue points towards A(1) or B(-1) and if the cue is valid,i.e. cue_validity v > 0.5
      ## if validity v < 0.5, the cue values are inverted
      if (testset[r,cue_order[v]]== 1 & cue_validities[cue_order[v]] > 0.5){           
        ttb_predict[r] <- 1  # A(1) is better than B(0)
        
      } else if (testset[r,cue_order[v]]== 1 & cue_validities[cue_order[v]] < 0.5){ 
        ttb_predict[r] <- -1  # B(-1) is better than A(1)
        
      } else if (testset[r,cue_order[v]]== -1 & cue_validities[cue_order[v]] > 0.5){
        ttb_predict[r] <- -1 # B(-1) is better than A(1)
        
      } else if (testset[r,cue_order[v]]== -1 & cue_validities[cue_order[v]] < 0.5){
        ttb_predict[r] <- 1  # A(1) is better than B(0)  
      }           
    } # end of trials for loop
    
    pred.ttb[[i]] <- ttb_predict # vector of length testset
    
  } # end of k for loop
  
  return(pred.ttb) # return the whole list with 20 elements
  
} # end of function




######  11.) TTB graph function (proportion correct) ######################


ttb.graph <- function(predictions.ttb, test.labels)
{
  
  ## function provides a summary of the performance (accuracy) of the tallying 
  ## heuristic
  ## input 'predictions.tbb' is always only 1 element of the list (vector)
  
  n <- length(predictions.ttb)
  m <- length(test.labels)
  
  prop.ttb <- sum(test.labels==predictions.ttb)/m
  
  # prop contains accuracies
  return(prop.ttb) # a single number
  
}