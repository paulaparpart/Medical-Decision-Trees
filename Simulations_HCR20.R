
####################################################################
#  R program that simulates people's decisions in a binary decision task
#  with multiple predictors and multiple environments. The environments 
#  vary in terms of: uncertainty (noise levels), sparsity structure
#	(of the regression weights), the number of predictors, sample size N, 
#	covariance among predictors,  
#	
#	The models used to predict future decisions are: 
# 1) Bayesian ridge regression (sampling from a truncated Gaussian posterior) is cross validated, but
# using only positive betas is equivalent to non-crossvalidated tallying heuristic
#	2) Tallying heuristic not cross-validated (assuming cue directionalities are known) 
# 3) Truncated Bayesian linear regression: sampling from Bayesian posterior with penalty=0 is linear
# regression model 
# 4) TTB cross-validated, assuming cue orders have to be learned from samples (Training sets) first
# 5) Lasso regression, assuming cue validities have to be learned from samples, i.e. no truncation, 
#   regression weights are derived from fitting regular lasso to training data

#	Model parameters: penalty parameter, N
#
#	Cross-Validation Simulation: Bayesian Ridge Regression vs. Heuristics
#     ========================================================
####### HCR-20 DATA SET Cross-validation  ####################################

## comparing logistic regression and logistic ridge regression 
## getting the optimal weights from best performing penalty parameter in the test set


############################################################################

#  	load libraries
require(clusterGeneration)
require(Matrix)
library(MASS)
require(ridge)
library(scatterplot3d)
library(msm)
library(tmvtnorm)
library(lars)
library(glmnet)
library(mcmc)
library(nnet)
library(xlsx)
library(e1071)


# source all functions including ridge
source("Functions_HCR20.R")

# load("HCR_results_new.RData") 
# load("HCR.eitheror.RData") 
# load("HCR.genderage.RData") 

start <- Sys.time()
seed <- 30

runme <- function(seed)
{
  
  ## This function is called when the user wants to generate some
  ## cross-validation analysis of artificial binary data fitting several differen
  ## models
  ##
  ## Arguments:
  ##
  ## seed 	A seed parameter has to be provided as an integer to the function.
  ##		The seed chosen for this simulation is 20. It is assigned at the 
  ##		bottom of the script when the runme function is called. The seed 
  ##		is used to generate the random choice of cross-validation 
  ##		partition. Using the same seed again will give the same results.
  

  ## 1.) READ-IN DATA #################################
  
  # 17 data sets overall available for simulation
#    for (e in 1:2){
    
    e <- 1       
    worlds <- c("physical_aggression", "verbal_aggression", "Either_or_aggression")
    env <- worlds[e]
       
    if (e==1)  {
      dataset <- read.table("physical_aggression_genderage.txt",header=TRUE)
    } else if (e==2)     {
      dataset <- read.table("verbal_aggression_genderage.txt",header=TRUE)
    }
#     } else if (e==3)     {
#       dataset <- read.table("Either_or_aggression.txt",header=TRUE)
#     } 
 
    # to include gender and age, use:
    # physical_aggression_genderage.txt
    # verbal_aggression_genderage.txt
  
    ## to exclude gender and age, use:
    # physical_aggression.txt
    # verbal_aggression.txt


     #load("HCR.svm_lassoVar.RData") 
# 
#        if (e > 1){
#          load("HCR.svm_Nogenderage.RData")    
#         } 

    ## CHECK whether DV is in column 4 for every dataset!
    y.pos <- 4 ## moved  one further along to include HCR_Total!!
    
    # attach to be able to access headers
    attach(dataset)
    labels <- names(dataset)
    if(labels[y.pos+1] != 'Completeness_of_data'){
      stop('the dataset is shifted, completeness check is not in column 4!')
    }
    
    # cue are still starting 2 behind y.pos and til the last column of dataset
    # col_cues <- c(6:ncol(dataset))
    

    ## LASSO -selected variables ONLY for Physical aggression: 
## with gender/age: H1, H3, H5, C2,C3,C4,R1,R2,R4, R5, Gender, Age 
    
## THE BEST 10 PREDICTORS: 
    col_cues <- c(6,8,10,16,17,18,20,21,23,24,25,26)

## make sure that Age/Gender is last 2 columns and not summed in HCR_total score!



## THE WORST 9 PREDICTORS: 
    #col_cues <- c(7,9,11,12,13,14,15,19,22)
 


    #grep("R1", colnames(dataset))

## without gender/age: H1, H2, H3, H5, H6, H9,H10,C2,C3,C4,R1,R2,R4,R5 (no gender and age).
    #col_cues <- c(6,7,8,10,11,13,14,16,17,18,20,21,23,24)


    ## assess the number of predictor X's, from column 5 to the end
    Predictors <- length(col_cues)
     
    # number of objects
    N <- nrow(dataset)
    
    ## Choose the number of partitions for cross-validation
    k <- 100   # to test small training sets
    
    # -----------------------------------------------------------------------------------------------------
    ## BEFORE: average redundancy between cues
    cov <- as.matrix(dataset[,col_cues])
    cor_mat <- cor(cov,method = "pearson")
    # the lower triangle of the cov matrix contains pairwise correlations
    min_cor <- min(cor_mat[lower.tri(cor_mat)])
    max_cor <- max(cor_mat[lower.tri(cor_mat)])
    abs_mean_cor <- mean(abs(cor_mat[lower.tri(cor_mat)]))
    # -------------------------------------------------------------------------------------------------  
    
    ###-------------------- Create the cue matrix with dependent in a data frame  ######################################
    ##no binary comparisons, just original cue data   
    bdata_diff <- dataset[ ,col_cues] # is data frame already
    dependent <- dataset[  ,y.pos]   
    paired_data <- data.frame(bdata_diff, dependent)
    paired_data$dependent <- as.factor(paired_data$dependent) # make sure y is a factor for svm
    ## -----------------------------------------------------------------------------------------
    
##--------------- get the logistic regression coefficients of the whole data set --------------------------------------------------
          d.pos <- ncol(paired_data)
          variables <- names(paired_data) # of training set only!
          fmla <- as.formula(paste(variables[d.pos]," ~ ", paste(variables[-d.pos], collapse= "+"), paste(" -1")))
          
          logistic.overall <- glm(fmla, family=binomial, data=paired_data)
          summary(logistic.overall)
          logistic_coefficients <- logistic.overall$coeff

# also print the signifiance p values if possible!
### -----------------------------------------------------------------------------------

      #### --------	Generate the 1st cross-validation partitions to estimate optimal parameters
      seed <- 30 # for permuatation
      #hold-one-out
      percent <- 0.10# test set size!  # ORDER: 0.10, 0.25,0.50,0.75, 0.80 
      percent_training <- (1 - percent)  # training is the complement!
      training_items <- (1-percent)*N*k    
      cv <- cv.indexing(k,nrow(paired_data),percent)
 
      ### --------------- Get the 2nd cross-validation sets to assess performance with optimal parameters
      seed <- 30 # for permuatation
      #hold-one-out
      percent <- 0.10 # test set size!
      percent_training <- (1 - percent)  # training is the complement!
      training_items <- (1-percent)*N*k
      # is different than above
      cv2 <- cv.indexing(k,nrow(paired_data),percent)
      

##------------------- Model SVM: Support vector Machine with Gaussian Kernel (RBF) ----------------------------------------###
      
      ## 2 parameters to tune
      gamma <- 2^(seq(-15,3,2))
      Cost <- 2^(seq(-5,15,2))
      
      ## Function cross-validating the svm, tuning 2 parameters on cv and predicting new data in cv2
      #y.pos <- ncol(paired_data)
      predictions.svm <- svm.method(paired_data, cv, cv2, ncol(paired_data), gamma, Cost) 

      # predictions.svm[[1]] contains a list with k elements, all k predictions vectors for optimal parameters
      svm_parameters <- predictions.svm[[2]] # contains optimal_parameters
      # predictions.svm[[3]] contains all support vector coefficients 

      svm_predictions <- predictions.svm[[1]]
      results.svm <- vector('numeric', k)
      for (i in 1:k) {
        testset2 <-  paired_data[cv2[i,]==0, ]
        test.labels <- testset2[ ,ncol(paired_data)]      
        results.svm[i] <- regression.graph(svm_predictions[[i]], test.labels) #equal to function results!
      }
      acc.svm <- sum(results.svm)/k


###--- ------------------  Model 0): TRADITIONAL LOGISTIC REGRESSION  ------------------------------------##################################### 

## dependent is binary variables: 1) Physical aggression (0/1) or Verbal aggression (0/1)

   #y.pos <- ncol(paired_data)
    predictions.regression <- logistic.regression(paired_data, cv, cv2, ncol(paired_data)) 
    
# predictions.regression[[1]] contains a list with k elements, all k predictions vectors for optimal threshold cut-off

    ####    Regression accuracy (generalization performance) ###########################
    results.regression <- vector('numeric', k)
    for (i in 1:k) {
      testset <-  paired_data[cv2[i,]==0, ] 
      test.labels <- testset[ ,ncol(paired_data)]      
      results.regression[i] <- regression.graph(predictions.regression[[i]], test.labels)
    }
    acc.regression <- sum(results.regression)/k

    ###--- ------------------  Model 1): TALLYING HEURISTIC ------------------------------------##################################### 
#     
#     predictions.tallyLearn <- tallying.learning(paired_data, cv2, ncol(paired_data) , Predictors)
#     
#     ####    Tallying accuracy (generalization performance) ###########################
#     results.tallying <- vector('numeric', k)
#     for (i in 1:k) {
#       testset <-  paired_data[cv[i,]==0, ]
#       test.labels <- testset[ ,ncol(paired_data)]      
#       results.tallying[i] <- tallying.graph(predictions.tallyLearn[[i]], test.labels)
#     }
#     acc.tallying <- sum(results.tallying)/k
#     
#     
#     ###--- ------------------  Model 2): Take-The-Best (TTB) HEURISTIC ------------------------------------#####################################
#     
#     predictions.ttb <- ttb.predictions(paired_data, cv2, ncol(paired_data), Predictors)
#     # predictions are all -1/+1 (like y!) 
#     
#     ####### TTB accuracy (generalization performance) ######################
#     results.ttb <- vector('numeric', k)
#     for (i in 1:k) {
#       testset <-  paired_data[cv[i,]==0, ]
#       test.labels <- testset[ ,ncol(paired_data)]
#       results.ttb[i] <- ttb.graph(predictions.ttb[[i]], test.labels)
#     }
#     # overall TTB accuracy across 1000 testsets
#     acc.ttb <- sum(results.ttb)/k  
#     
    #--------------------------------------------------------------------------------------------------------   
    
    
    ###--- ------------------  Model 3): Regular Ridge Regression (L2) ------------------------------------#####################################
    # set penalty paramter range
    penalty <- sort(c(0,exp(seq(0:30)-11),700))  # 32 penalty paramters plus a few more before 1000
    #penalty <- penalty[c(1,10,12,14,15,17,18,20,22,25,27,30,33)] # shorter version
    
    predictions.regridge <- regular.ridge(paired_data, cv, cv2, ncol(paired_data), penalty) 
    
    # predictions.regridge[[1]] contains a list with k elements, all k predictions vectors for optimal penalty
    # predictions.regridge[[2]] contains optimal_p
    # predictions.regridge[[3]] contains optimal coefficients optimal_coeff
      predictions.ridge <- predictions.regridge[[1]]
      results.ridge <- vector('numeric', k)
      for (i in 1:k) {
        testset <-  paired_data[cv2[i,]==0, ]
        test.labels <- testset[ ,ncol(paired_data)]      
        results.ridge[i] <- regression.graph(predictions.ridge[[i]], test.labels)
      }
      acc.ridge <- sum(results.ridge)/k



    ###--- ------------------  Model 4): Regular Lasso Regression (L1) ------------------------------------####################################
    penalty <- sort(c(0,exp(seq(0:30)-11),700))
    
    predictions.lasso <- regular.lasso(paired_data, cv, cv2, ncol(paired_data), penalty) 
    
    # predictions.lasso[[1]] contains a list with k elements, all k predictions vectors for optimal penalty
    # predictions.lasso[[2]] contains optimal_p
    # predictions.lasso[[3]] contains optimal coefficients optimal_coeff  
    predictions.l <- predictions.lasso[[1]]
    results.lasso <- vector('numeric', k)
    for (i in 1:k) {
      testset <-  paired_data[cv2[i,]==0, ]
      test.labels <- testset[ ,ncol(paired_data)]      
      results.lasso[i] <- regression.graph(predictions.l[[i]], test.labels)
    }
    acc.lasso <- sum(results.lasso)/k



    ###--- ------------------  Model 5): Simple Scale (optimal Threshold) ------------------------------------####################################
      
      ## Create data frame for simple scale model
      ## Total score is always in column 3  
    



        hcr_total <- rowSums(dataset[ ,col_cues[1:(length(col_cues)-2)]])

        dependent <- dataset[ ,4] 
        scale_data <- data.frame(hcr_total = hcr_total, dependent)
          
        predictions.scale <- simple.scale(scale_data, cv, cv2, ncol(scale_data), Predictors) 
      # returns a list predictions.scale[[i]] containing predictions vector
      
        ####   Simple scale accuracy (generalization performance) ###########################
         
           predictions.only  <- predictions.scale[[1]]
          results.scale <- vector('numeric', k)
          for (i in 1:k) {            
            testset <- scale_data[cv2[i,]==0, ]    
            test.labels <- testset[ ,2]  
            results.scale[i] <- regression.graph(predictions.only[[i]], test.labels)
          }
          acc.scale <- sum(results.scale)/k
         

   ##------------  Accordance between Logistic Regression and LASSO\Ridge ----------------------------------------------
    
    #acc_logisticRidge <- vector("numeric",k)
    #acc_logisticLasso <- vector("numeric",k)
    
    ## Agreement of Scale with new methods!!!
    acc_logisticScale <- vector("numeric",k)
    acc_ScaleRIdge  <- vector("numeric",k)
    acc_ScaleLasso  <- vector("numeric",k)

    ridge_predictions <-  predictions.regridge[[1]] # contains kxpredictions matrix of ridge
    ## ridge_matrix[i, ] # is numeric predictions vector per test set i for comparison    
    lasso_predictions <- predictions.lasso[[1]]
    ## lasso_matrix[i, ] # is predictions vector per test set i for comparison
    
    #regression.predictions is a list/ only necessary with optimal threshold!
#     regression.predictions <- predictions.regression[[1]]
    predictions.only  <- predictions.scale[[1]]

    for (i in 1:k){ 

      #acc_logisticRidge[i] <- accordance(predictions.regression[[i]],ridge_matrix[i, ])
      #acc_logisticLasso[i] <- accordance(predictions.regression[[i]],lasso_matrix[i, ])   
      acc_logisticScale[i] <- accordance(predictions.regression[[i]],predictions.only [[i]])
      acc_ScaleRIdge[i]  <- accordance(ridge_predictions[[i]],predictions.only [[i]])
      acc_ScaleLasso[i] <- accordance(lasso_predictions[[i]],predictions.only [[i]])
      
    }
    

    #Log_Ridge <- sum(acc_logisticRidge)/k
    #Log_Lasso <- sum(acc_logisticLasso)/k
    Log_Scale <- sum(acc_logisticScale)/k
    Scale_Ridge <-sum(acc_ScaleRIdge)/k
    Scale_Lasso <-sum(acc_ScaleLasso)/k

    ###################### STORING RESULTS ############################################
    
    ## the standard deviation of the average of 1000 runs, if low then good ??
    ridge_coeff <- as.matrix(predictions.regridge[[3]])
    lasso_coeff <- as.matrix(predictions.lasso[[3]])                           
    
# predictions.lasso[[1]] contains a list with k elements, all k predictions vectors for optimal penalty
# predictions.lasso[[2]] contains optimal_p
# predictions.lasso[[3]] contains optimal coefficients optimal_coeff  

    # put all outcome variables into a data frame
    hcr.svm <- data.frame(environment = env,
                              acc_scale = round(acc.scale,3), acc_regression = round(acc.regression,3), 
                              acc_ridge = round(acc.ridge,3),acc_lasso = round(acc.lasso,3),
                              acc_svm = round(acc.svm,3),
                              penalty_ridge = predictions.regridge[[2]],         
                              penalty_lasso = predictions.lasso[[2]], param_svm.Cost =svm_parameters$Cost,
                              param_svm.gamma =svm_parameters$gamma,
                              Scale_Regression =  Log_Scale, Scale_Ridge = Scale_Ridge,Scale_Lasso = Scale_Lasso,
                              Predictors, N,k,percent_training = percent_training,
                              training_items = training_items, cov_mean = abs_mean_cor,
                              cov_min = min_cor, cov_max = max_cor,  
                              logistic_coefficients = logistic_coefficients,
                              ridge_regularized = ridge_coeff[ ,1],lasso_regularized =  lasso_coeff[ ,1])
    

    print(hcr.svm)

    HCR.svm_lassoVar_Best <- hcr.svm


    #HCR.svm_lassoVar <- rbind(HCR.svm_lassoVar, hcr.svm) 
   
#     # combine new and old data (old data first, then new data below) only if e==2
#     if (e == 1){
#        HCR.svm_Nogenderage <- hcr.svm  
#     } else {
#       HCR.svm_Nogenderage <- rbind(HCR.svm_Nogenderage, hcr.svm) 
#     } 
    
    # save R data frame
    save(HCR.svm_lassoVar_Best,file="HCR.svm_lassoVar_Best.RData") 
    # write excel workbook 
    write.xlsx(HCR.svm_lassoVar_Best, file="HCR.svm_lassoVar_Best.xlsx", row.names=TRUE)
    
    detach(dataset)

    
#   } # environments loop 
# 

}


results <- runme(seed)


end <- Sys.time()
duration <- end - start
cat(paste("\n\nSIMULATION DURATION  = ",duration))

