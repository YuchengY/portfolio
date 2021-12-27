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


# register cores, initialize clusters and set seed

ncores <- 12

set.seed(4352216)

#just two sub folders this time? 
#"Result/Special" and "Result/Lin_nonnormal"

folders <- paste0("Results2/", dir("Results2"))

analysis_results <- map_dfr(folders, ~get_test_results(.x, ncores = ncores))

write.csv(analysis_results, "test_results.csv")
