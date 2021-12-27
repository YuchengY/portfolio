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
# library(doParallel)
# library(foreach)
library(parallel) ## < for mclapply
library(purrr)
library(readr)


# Call Functions

source("Helper_Functions.R")


# register cores, initialize clusters and set seed

ncores <- ncores <- detectCores() - 1

set.seed(4352216)


folders <- paste0("Result/", dir("Result"))

analysis_results <- map_dfr(folders, ~ get_test_results(.x, ncores = ncores))

write.csv(analysis_results, "test_results.csv")
