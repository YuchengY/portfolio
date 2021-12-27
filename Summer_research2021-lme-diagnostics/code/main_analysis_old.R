#!/usr/bin/env /opt/R/4.0.4/bin/Rscript

# This script runs tests on residuals from simulated models

# load packages
library(lmtest)
library(dplyr)
library(tidyr)
library(lme4)
library(nlme)
library(HLMdiag)
library(PearsonDS)
library(MASS)
#library(doParallel)
#library(foreach)
library(parallel) ##< for mclapply
library(purrr)
library(readr)


# Call Functions

source("Helper_Functions.R")
source("Simulation_Functions.R")


# register cores, initialize clusters and set seed

ncores <- 12      # 12 core / Bio 328
cl <- makeCluster(ncores, type = "FORK")
# registerDoParallel(cl)

set.seed(4352216)
nsim <- 1000        # nsim Need to be consistent with main_resid.R !! (change to 5 if doing test)



# Analysis and Testing 


# Important Helpers (Already Put in Helper_Functions.R)

# calculate rate of p < 0.05
rate <- function (df){
  mean(df < 0.05)
}

# conduct shapiro, bp, anova, and normality test of re
all_residual_tests<- function(model_residual, design_matrix){
  #norm 
  tests <- 
    foreach(j = 1: nrow(design_matrix))%:%
    foreach(i = 1:nsim, .combine=c ) %dopar% {   # add .combine = c
      setting <- model_residual[[j]]
      
      #norm
      shapiro<-shapiro.test(setting[[i]]$first_resid$.std.resid)$p.value
      
      #hetero
      bp_test<-bp(setting[[i]]$first_resid$.resid, setting[[i]]$first_resid$.fitted)["p_val"]
      
      #re_norm
      re_norm<-ks.test(setting[[i]]$second_resid$mahalanobis, "pchisq", df=2)$p.value
      
      #linearity 
      chol.mar.resid <- setting[[i]]$first_resid$.chol.mar.resid
      fittted  <-  setting[[i]]$first_resid$.mar.fitted
      lin_test <- anova(lm(chol.mar.resid ~ fittted),
                        lm(chol.mar.resid ~ fittted + I(fittted^2)))$`Pr(>F)`[2]
      
      #bind them together and give names
      all_tests <- list(shapiro, bp_test, re_norm, lin_test)
      names(all_tests) <- c("shapiro", "bp_test", "re_norm", "lin_test")
      all_tests
    }
  
  nms <-c("shapiro", "bp_test.p_val", "re_norm", "lin_test")
  
  reject_rate <- 
    foreach( j= 1:nrow(design_matrix)) %:%
    foreach (i = 1:length(nms), .combine = c) %dopar%{
      mean(unlist(tests[[j]])[names(unlist(tests[[j]])) == as.character(nms[i])] < 0.05)
    }
  
  return (reject_rate)
}


#test results df ready for adding into design matrix 
test.df <- function(test_results, design_matrix){
  nms <-c("shapiro", "bp_test.p_val", "re_norm", "lin_test")
  df<- foreach (i = 1:length(nms)) %:%
    foreach (j = 1:nrow(design_matrix))%dopar% {
      test_results[[j]][i]
    }
  return (df)
}


# conduct four test and save with design matrix
conduct_4_test <- function(model_residuals, design_matrix){
  
  test_results <- all_residual_tests(model_residuals, design_matrix)
  
  df <- test.df(test_results, design_matrix)
  
  design_results<-design_matrix
  design_results$shapiro<-as.numeric(df[[1]])
  design_results$bp_test<-as.numeric(df[[2]])
  design_results$re_norm<-as.numeric(df[[3]])
  design_results$lin_test<-as.numeric(df[[4]])
  
  return(design_results)
}




###-----------------------------------

#now it's the mclapply version of the above functions

all_residual_tests_mc<- function(model_residual, design_matrix){
  #norm 
  tests <- 
    mclapply(1: nrow(design_matrix), function(j){ 
      mclapply(1:nsim, function (i){
        setting <- model_residual[[j]]
        
        #norm
      shapiro<-shapiro.test(setting[[i]]$first_resid$.std.resid)$p.value
      
      #hetero
      bp_test<-bp(setting[[i]]$first_resid$.resid, setting[[i]]$first_resid$.fitted)["p_val"]
      
      #re_norm
      re_norm<-ks.test(setting[[i]]$second_resid$mahalanobis, "pchisq", df=2)$p.value
      
      #linearity 
      chol.mar.resid <- setting[[i]]$first_resid$.chol.mar.resid
      fittted  <-  setting[[i]]$first_resid$.mar.fitted
      lin_test <- anova(lm(chol.mar.resid ~ fittted),
                        lm(chol.mar.resid ~ fittted + I(fittted^2)))$`Pr(>F)`[2]
      
      #bind them together and give names
      all_tests <- list(shapiro, bp_test, re_norm, lin_test)
      names(all_tests) <- c("shapiro", "bp_test", "re_norm", "lin_test")
      all_tests
      }, mc.cores = ncores)}, mc.cores = ncores
      )
  
    nms <-c("shapiro", "bp_test.p_val", "re_norm", "lin_test")

  reject_rate <-
    mclapply(1: nrow(design_matrix), function(j){
      mclapply(1:length(nms), function (i){
      mean(unlist(tests[[j]])[names(unlist(tests[[j]])) == as.character(nms[i])] < 0.05)
    }, mc.cores =ncores)}, mc.cores =ncores)

 
  return (reject_rate)
}

#test.df stays the same


conduct_4_test_mc <- function(model_residuals, design_matrix){
  
  test_results <- all_residual_tests_mc(model_residuals, design_matrix)
  
  df <- test.df(test_results, design_matrix)
  
  design_results<-design_matrix
  design_results$shapiro<-as.numeric(flatten(df[[1]]))
  design_results$bp_test<-as.numeric(flatten(df[[2]]))
  design_results$re_norm<-as.numeric(flatten(df[[3]]))
  design_results$lin_test<-as.numeric(flatten(df[[4]]))
  
  return(design_results)
}



# Need every design matrix 
# Currently 108 rows 

variance<-c("large_error","large_re") # var settings added 7.25
balance<-c("same","balanced","unbalanced")
norm<-c("norm", "skewness_3", "skewness_1.5", "skewness_0.8", "bimodal")
hetero<-c("0","2","4","8")
lin<-c("linear","sq")
omit<-c("full","reduced")
RE <- c("norm_re","mildly_skewed_re_intercept", "mildly_skewed_re_slope")
hetero_lin <- c("linear_homo","sq_2","sq_4","sq_8")
hetero_norm <- c("0_skew","2_skew", "4_skew", "8_skew")

# one good and seven other scenarios
design_matrix_good<-expand.grid(variance, balance, norm[1], hetero[1],lin[1],RE[1],omit[1],hetero_lin[1],hetero_norm[1])
design_matrix_norm<-expand.grid(variance, balance, norm[-1], hetero[1],lin[1],RE[1],omit[1],hetero_lin[1],hetero_norm[1])
design_matrix_heter<-expand.grid(variance, balance, norm[1], hetero[-1],lin[1],RE[1],omit[1],hetero_lin[1],hetero_norm[1])
design_matrix_lin<-expand.grid(variance, balance, norm[1], hetero[1],lin[-1],RE[1],omit[1],hetero_lin[1],hetero_norm[1])
design_matrix_re<-expand.grid(variance, balance, norm[1], hetero[1],lin[1],RE[-1],omit[1],hetero_lin[1],hetero_norm[1])
design_matrix_omit<-expand.grid(variance, balance, norm[1], hetero[1],lin[1],RE[1],omit[-1],hetero_lin[1],hetero_norm[1])
design_matrix_hetero_lin<-expand.grid(variance, balance, norm[1], hetero[1],lin[1],RE[1],omit[1],hetero_lin[-1],hetero_norm[1])
design_matrix_hetero_norm<-expand.grid(variance, balance, norm[1], hetero[1],lin[1],RE[1],omit[1],hetero_lin[1],hetero_norm[-1])

# Rename colnames
scenario <- c("variance","balance", "normality", "heteroscedasticity", "linearity", "random_effect", "fixed_effect","Hetero_lin","Hetero_norm")
colnames(design_matrix_norm)<- colnames(design_matrix_heter) <- 
  colnames(design_matrix_lin) <- colnames(design_matrix_re)  <-
  colnames(design_matrix_omit) <- colnames(design_matrix_hetero_lin) <-
  colnames(design_matrix_hetero_norm) <- colnames(design_matrix_good) <- scenario

# Combine all settings considered
design_matrix<- rbind(design_matrix_good,
                      design_matrix_norm,
                      design_matrix_heter,
                      design_matrix_lin,
                      design_matrix_re,
                      design_matrix_omit,
                      design_matrix_hetero_lin,
                      design_matrix_hetero_norm) %>% distinct()






## Scenario 0 Residuals from Good reference models

good_residuals <- 
  mclapply( 1:nrow(design_matrix_good), function (i){
    readRDS(paste("Result/Good/scenario_0_", as.character(i), ".RDS", sep = ""))
  },  mc.cores = ncores
)


design_matrix_good_results <- conduct_4_test_mc(good_residuals, design_matrix_good)




## Scenario 1 Normaltiy test

#read all RDS as a nested list
norm_residuals <-
  mclapply( 1:nrow(design_matrix_norm), function (i){ readRDS(paste("Result/Normality/scenario_1_", as.character(i), ".RDS", sep = ""))
  }, mc.cores = ncores
)


#1. conduct shapiro test

# Norm_residual_test <- function(model_residual){
#   
#   norm_tests <- 
#     foreach(j = 1: nrow(design_matrix_norm))%:%
#     foreach(i = 1:nsim, .combine = c) %dopar% {   # add .combine = c
#     setting <- model_residual[[j]]
#     shapiro.test(setting[[i]]$first_resid$.std.resid)$p.value
#     }
#     
#   reject_rate<-lapply(norm_tests, rate)
#  
#  
#   return (reject_rate)
# }

#rejection rate attached to design matrix
# norm_reject <- Norm_residual_test(norm_residuals)
# shapiro_norm_results<-design_matrix_norm
# shapiro_norm_results$results<-as.numeric(norm_reject)


# 2. Four test of this scenario
design_matrix_norm_results <- conduct_4_test_mc(norm_residuals, design_matrix_norm)



## Scenario 2 Heteroscedasticity test

## Read residuals for constant variance analysis

hetero_residuals <-
  mclapply( 1:nrow(design_matrix_heter), function (i){
    readRDS(paste("Result/Heteroscedasticity/scenario_2_", as.character(i), ".RDS", sep = ""))
  }, mc.cores = ncores
)



design_matrix_hetero_results <- conduct_4_test_mc(hetero_residuals, design_matrix_heter)



## Scenario 3 Residual Linearity

## Read residuals for linearity analysis
linearity_residuals <-
   mclapply( 1:nrow(design_matrix_lin), function (i){
    readRDS(paste("Result/Linearity/scenario_3_", as.character(i), ".RDS", sep = ""))
   }, mc.cores = ncores
)


design_matrix_linear_results <- conduct_4_test_mc(linearity_residuals, design_matrix_lin)



## Scenario 4 RE

#read RDS as a list
re_residuals <-
  mclapply( 1:nrow(design_matrix_re), function (i){
    readRDS(paste("Result/Randomeffect/", "scenario_4_", as.character(i), ".RDS", sep = ""))
  }, mc.cores = ncores
)



design_matrix_re_results <- conduct_4_test_mc(re_residuals, design_matrix_re)



## Scenario 5 Hetero + Nonnormal  

#extract residuals 
hetero_norm_residuals <-
  mclapply( 1:nrow(design_matrix_hetero_norm), function (i){
    readRDS(paste("Result/Heter_nonnormal/", "scenario_5_", as.character(i), ".RDS", sep = ""))
  }, mc.cores = ncores
)


#conduct bp test 

#helper
# hetero_norm_residual_test <- function(model_residual){
#   
#   tests <- 
#     foreach(j = 1: nrow(design_matrix_hetero_norm))%:%
#     foreach(i = 1:nsim, .combine = c) %dopar% {
#     setting <- model_residual[[j]]
#     bp(setting[[i]]$first_resid$.resid, setting[[i]]$first_resid$.fitted)["p_val"]
# 
#     }
# 
#   reject_rate<-lapply(tests, rate)
#    
#   return (reject_rate)
# }


#rejection rate attached to design matrix

# hetero_norm_reject<-hetero_norm_residual_test(hetero_norm_residuals)
# hetero_norm_test_results<-design_matrix_hetero_norm
# hetero_norm_test_results$results<-as.numeric(hetero_norm_reject)


# Four test on this scenario
design_matrix_hetero_norm_results <- conduct_4_test_mc(hetero_norm_residuals, design_matrix_hetero_norm)




## Scenario 6 Hetero + Nonlinaer

hetero_lin_residuals <-
  mclapply( 1:nrow(design_matrix_hetero_lin), function (i){
    readRDS(paste("Result/Heter_nonlinear/", "scenario_6_", as.character(i), ".RDS", sep = ""))
  }, mc.cores = ncores
)


design_matrix_hetero_lin_results <- conduct_4_test_mc(hetero_lin_residuals, design_matrix_hetero_lin)



## Scenario 7 Omitting fixed effect

omit_residuals <-
  mclapply( 1:nrow(design_matrix_omit), function (i){
    readRDS(paste("Result/Fixedeffect/", "scenario_7_", as.character(i), ".RDS", sep = ""))
  }, mc.cores = ncores
)


design_matrix_omit_results <- conduct_4_test_mc(omit_residuals, design_matrix_omit)




# All in one 
result_matrix<- rbind(design_matrix_good_results,
                      design_matrix_norm_results,
                      design_matrix_hetero_results,
                      design_matrix_linear_results,
                      design_matrix_re_results,
                      design_matrix_hetero_norm_results,
                      design_matrix_hetero_lin_results,
                      design_matrix_omit_results) %>% distinct()

result_matrix <- as_tibble(result_matrix)

# Save as csv in the folder
write.table(result_matrix, file = "test_results_old.csv", sep = ",", col.names = NA,
            qmethod = "double")




# stop cluster
stopCluster(cl)

