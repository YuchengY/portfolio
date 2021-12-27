#!/usr/bin/env /opt/R/4.0.4/bin/Rscript

# This script extracts residuals from simulated models

# load packages
library(lmtest)
library(dplyr)
library(tidyr)
library(lme4)
library(nlme)
library(HLMdiag)
library(PearsonDS)
library(MASS)
library(parallel)
library(purrr)
library(readr)

# Call Functions

source("Helper_Functions.R")
source("Simulation_Functions.R")


# register cores, initialize clusters and set seed

ncores <- detectCores() - 1 # Use 12 core / Bio 328
cl <- makeCluster(ncores, type = "FORK")

set.seed(4352216)
nsim <- 1000 # change as needed


# Design matrix

variance <- c("large_error", "large_re")
balance <- c("same", "balanced", "unbalanced")
norm <- c("norm", "skewness_3", "skewness_1.5", "skewness_0.8", "bimodal")
hetero <- c("0", "2", "4", "8")
lin <- c("linear", "sq")
omit <- c("full", "reduced")
RE <- c("norm_re", "mildly_skewed_re_intercept", "mildly_skewed_re_slope")
hetero_lin <- c("linear_homo", "sq_2", "sq_4", "sq_8")
hetero_norm <- c("0_skew", "2_skew", "4_skew", "8_skew")
lin_norm <- c("linear_norm", "reduced_skew", "reduced_bimodal") # mildly skewed, skewness = 1.5
special <- c("standard", "time_seq", "ar_error", "re_int")

# one good and seven other scenarios
design_matrix_good <- expand.grid(variance, balance, norm[1], hetero[1], lin[1], RE[1], omit[1], hetero_lin[1], hetero_norm[1], lin_norm[1], special[1])
design_matrix_norm <- expand.grid(variance, balance, norm[-1], hetero[1], lin[1], RE[1], omit[1], hetero_lin[1], hetero_norm[1], lin_norm[1], special[1])
design_matrix_heter <- expand.grid(variance, balance, norm[1], hetero[-1], lin[1], RE[1], omit[1], hetero_lin[1], hetero_norm[1], lin_norm[1], special[1])
design_matrix_lin <- expand.grid(variance, balance, norm[1], hetero[1], lin[-1], RE[1], omit[1], hetero_lin[1], hetero_norm[1], lin_norm[1], special[1])
design_matrix_re <- expand.grid(variance, balance, norm[1], hetero[1], lin[1], RE[-1], omit[1], hetero_lin[1], hetero_norm[1], lin_norm[1], special[1])
design_matrix_omit <- expand.grid(variance, balance, norm[1], hetero[1], lin[1], RE[1], omit[-1], hetero_lin[1], hetero_norm[1], lin_norm[1], special[1])
design_matrix_hetero_lin <- expand.grid(variance, balance, norm[1], hetero[1], lin[1], RE[1], omit[1], hetero_lin[-1], hetero_norm[1], lin_norm[1], special[1])
design_matrix_hetero_norm <- expand.grid(variance, balance, norm[1], hetero[1], lin[1], RE[1], omit[1], hetero_lin[1], hetero_norm[-1], lin_norm[1], special[1])
design_matrix_lin_norm <- expand.grid(variance, balance, norm[1], hetero[1], lin[1], RE[1], omit[1], hetero_lin[1], hetero_norm[1], lin_norm[-1], special[1])
design_matrix_spec <- expand.grid(variance, balance, norm[1], hetero[1], lin[1], RE[1], omit[1], hetero_lin[1], hetero_norm[1], lin_norm[1], special[-1])

# Rename colnames
scenario <- c("variance", "balance", "normality", "heteroscedasticity", "linearity", "random_effect", "fixed_effect", "Hetero_lin", "Hetero_norm", "lin_norm", "special")
colnames(design_matrix_norm) <- colnames(design_matrix_heter) <-
  colnames(design_matrix_lin) <- colnames(design_matrix_re) <-
  colnames(design_matrix_omit) <- colnames(design_matrix_hetero_lin) <-
  colnames(design_matrix_hetero_norm) <- colnames(design_matrix_spec) <-
  colnames(design_matrix_lin_norm) <- colnames(design_matrix_good) <- scenario

# Combine all settings considered
design_matrix <- rbind(
  design_matrix_good,
  design_matrix_norm,
  design_matrix_heter,
  design_matrix_lin,
  design_matrix_re,
  design_matrix_omit,
  design_matrix_hetero_lin,
  design_matrix_hetero_norm,
  design_matrix_lin_norm,
  design_matrix_spec
) %>% distinct()



# Get simualted Data list and Extract residuals




## Scenario 7
### Analysis of Model Residuals from models with/without omitting fix effect

fixed_data_list <- flattenlist(get_design_matrix_data(design_matrix_omit))

sim.fm.fix <- mclapply(
  fixed_data_list,
  FUN = function(x) {
    sim.hlmfix(x, z.dsn = as.character(x$omit[1]), var = as.character(x$variance[1]))
  }, mc.cores = ncores
)

fix_omit_resid <- extract_2levels_resid_omit(sim.fm.fix, design_matrix_omit, nsim)

for (k in 1:length(fix_omit_resid)) {
  saveRDS(fix_omit_resid[[k]], paste("Result3/Fixedeffect/", "scenario_7_", as.character(k), ".RDS", sep = ""))
}

rm(fixed_data_list)
rm(fix_omit_resid)




## Scenario 9
# Special Cases: AR, RE-intercept, Time-seq
spec_data_list <- unlist(get_design_matrix_data_2(design_matrix_spec), recursive = F)

sim.fm.spec <- mclapply(
  spec_data_list,
  FUN = function(x) {
    sim.hlm_spec(x$simdata.df, var = as.character(x$simdata.df$variance[1]), special = as.character(x$simdata.df$special[1]), J = x$J)
  }, mc.cores = ncores
)


# get resid
spec_resid <- extract_2levels_resid(sim.fm.spec, design_matrix_spec, nsim)

for (k in 1:length(spec_resid)) {
  saveRDS(spec_resid[[k]], paste("Result3/Special/", "scenario_9_", as.character(k), ".RDS", sep = ""))
}

rm(spec_data_list)
rm(spec_resid)




# Analysis & Test

folders <- paste0("Result3/", dir("Result3"))

analysis_results <- map_dfr(folders, ~ get_test_results_new(.x, ncores = ncores))

write.csv(analysis_results, "test_results_3.csv")



# stop cluster
stopCluster(cl)
