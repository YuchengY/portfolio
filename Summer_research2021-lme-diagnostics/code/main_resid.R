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
# library(doParallel)
# library(foreach)
library(parallel) ## < for mclapply
library(purrr)
library(readr)

# Call Functions

source("Helper_Functions.R")
source("Simulation_Functions.R")


# register cores, initialize clusters and set seed

ncores <- detectCores() - 1 # Use 12 core / Bio 328
cl <- makeCluster(ncores, type = "FORK")
# registerDoParallel(cl)

set.seed(4352216)
nsim <- 1000 # change as needed


# Design matrix

# New Idea about design matrix (Yicheng 7.24)
# Only include settings that we want to analyze in this study
# Currently 108 rows
variance <- c("large_error", "large_re") # var settings added 7.25
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
design_matrix_spec <- expand.grid(variance, norm[1], hetero[1], lin[1], RE[1], omit[1], hetero_lin[1], hetero_norm[1], balance, lin_norm[-1], special[1])
design_matrix_lin_norm <- expand.grid(variance, balance, norm[1], hetero[1], lin[1], RE[1], omit[1], hetero_lin[1], hetero_norm[1], lin_norm[1], special[-1])

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

## Scenario 0
### Analysis of Model Residuals from good models without any deliberate violation

good_data_list <- flattenlist(get_design_matrix_data(design_matrix_good))
# sim.fm.good<-foreach(i = 1:length(good_data_list)) %dopar%{
#   sim.hlm(good_data_list[[i]], e.dsn  = as.character(good_data_list[[i]]$normality[1]),
#           var = as.character(good_data_list[[i]]$variance[1]))
# }

## Not worth parallelizing here
sim.fm.good <- mclapply(
  good_data_list,
  FUN = function(x) {
    sim.hlm(x, e.dsn = as.character(x$normality[1]), var = as.character(x$variance[1]))
  }, mc.cores = ncores
)


reference_resid <- extract_2levels_resid(sim.fm.good, design_matrix_good, nsim)

for (k in 1:length(reference_resid)) {
  saveRDS(reference_resid[[k]], paste("Result/Good/", "scenario_0_", as.character(k), ".RDS", sep = ""))
}

rm(good_data_list)
rm(reference_resid)



## Scenario 1
### Analysis of Model Residuals from models with/without Normality violation

norm_data_list <- flattenlist(get_design_matrix_data(design_matrix_norm))

sim.fm.norm <- mclapply(
  norm_data_list,
  FUN = function(x) {
    sim.hlm(x, e.dsn = as.character(x$normality[1]), var = as.character(x$variance[1]))
  }, mc.cores = ncores
)


norm_resid <- extract_2levels_resid(sim.fm.norm, design_matrix_norm, nsim)

for (k in 1:length(norm_resid)) {
  saveRDS(norm_resid[[k]], paste("Result/Normality/", "scenario_1_", as.character(k), ".RDS", sep = ""))
}

rm(norm_data_list)
rm(norm_resid)



## Scenario 2
### Analysis of Model Residuals from models with/without Homoscedascity violation

hetero_data_list <- flattenlist(get_design_matrix_data(design_matrix_heter))

sim.fm.hetero <- mclapply(
  hetero_data_list,
  FUN = function(x) {
    sim.hlmvar(x, evar.dsn = as.numeric(x$hetero[1]), var = as.character(x$variance[1]))
  }, mc.cores = ncores
)

hetero_resid <- extract_2levels_resid(sim.fm.hetero, design_matrix_heter, nsim)

for (k in 1:length(hetero_resid)) {
  saveRDS(hetero_resid[[k]], paste("Result/Heteroscedasticity/", "scenario_2_", as.character(k), ".RDS", sep = ""))
}


rm(hetero_data_list)
rm(hetero_resid)



## Scenario 3
### Analysis of Model Residuals from models with/without Linearity violation

linear_data_list <- flattenlist(get_design_matrix_data(design_matrix_lin))

sim.fm.lin <- mclapply(
  linear_data_list,
  FUN = function(x) {
    sim.hlmlin(x, rel.dsn = as.character(x$linearity[1]), var = as.character(x$variance[1]))
  }, mc.cores = ncores
)

linear_resid <- extract_2levels_resid(sim.fm.lin, design_matrix_lin, nsim)

for (k in 1:length(linear_resid)) {
  saveRDS(linear_resid[[k]], paste("Result/Linearity/", "scenario_3_", as.character(k), ".RDS", sep = ""))
}


rm(linear_data_list)
rm(linear_resid)



## Scenario 4
### Analysis of Random effect from models with/without normality violation

random_effects_data_list <- flattenlist(get_design_matrix_data(design_matrix_re))

sim.fm.re <- mclapply(
  random_effects_data_list,
  FUN = function(x) {
    sim.hlm_re(x, re.dsn = as.character(x$random_effects[1]), var = as.character(x$variance[1]))
  }, mc.cores = ncores
)

RE_resid <- extract_2levels_resid(sim.fm.re, design_matrix_re, nsim)

for (k in 1:length(RE_resid)) {
  saveRDS(RE_resid[[k]], paste("Result/Randomeffect/", "scenario_4_", as.character(k), ".RDS", sep = ""))
}


rm(random_effects_data_list)
rm(RE_resid)



## Scenario 5
### Analysis of Model Residuals from models with/without Homoscedascity & normality

hetero_norm_data_list <- flattenlist(get_design_matrix_data(design_matrix_hetero_norm))

sim.fm.heteronorm <- mclapply(
  hetero_norm_data_list,
  FUN = function(x) {
    sim.hlmvar_norm(x, evar.dsn = as.character(x$hetero_norm[1]), var = as.character(x$variance[1]))
  }, mc.cores = ncores
)

hetero_norm_resid <- extract_2levels_resid(sim.fm.heteronorm, design_matrix_hetero_norm, nsim)

for (k in 1:length(hetero_norm_resid)) {
  saveRDS(hetero_norm_resid[[k]], paste("Result/Heter_nonnormal/", "scenario_5_", as.character(k), ".RDS", sep = ""))
}


rm(hetero_norm_data_list)
rm(hetero_norm_resid)



## Scenario 6
### Analysis of Model Residuals from models with/without Homoscedascity & Linearity

hetero_lin_data_list <- flattenlist(get_design_matrix_data(design_matrix_hetero_lin))

sim.fm.heterolin <- mclapply(
  hetero_lin_data_list,
  FUN = function(x) {
    sim.hlm_heterlin(x, helin.dsn = as.character(x$hetero_lin[1]), var = as.character(x$variance[1]))
  }, mc.cores = ncores
)

heterolin_resid <- extract_2levels_resid(sim.fm.heterolin, design_matrix_hetero_lin, nsim)

for (k in 1:length(heterolin_resid)) {
  saveRDS(heterolin_resid[[k]], paste("Result/Heter_nonlinear/", "scenario_6_", as.character(k), ".RDS", sep = ""))
}


rm(hetero_lin_data_list)
rm(heterolin_resid)



## Scenario 7
### Analysis of Model Residuals from models with/without omitting fix effect

fixed_data_list <- flattenlist(get_design_matrix_data(design_matrix_omit))

sim.fm.fix <- mclapply(
  fixed_data_list,
  FUN = function(x) {
    sim.hlmfix(x, z.dsn = as.character(x$omit[1]), var = as.character(x$variance[1]))
  }, mc.cores = ncores
)

fix_omit_resid <- extract_2levels_resid(sim.fm.fix, design_matrix_omit, nsim)

for (k in 1:length(fix_omit_resid)) {
  saveRDS(fix_omit_resid[[k]], paste("Result/Fixedeffect/", "scenario_7_", as.character(k), ".RDS", sep = ""))
}

rm(fixed_data_list)
rm(fix_omit_resid)



## Scenario 8
### Analysis of Model Residuals from models with/without Normality & Linearity
linnorm_data_list <- unlist(get_design_matrix_data_2(design_matrix_lin_norm), recursive = F)

sim.fm.linnorm <- mclapply(
  linnorm_data_list,
  FUN = function(x) {
    sim.hlm_linnorm(x$simdata.df, lin_norm.dsn = as.character(x$simdata.df$lin_norm[1]), var = as.character(x$simdata.df$variance[1]))
  }, mc.cores = ncores
)


lin_norm_resid <- extract_2levels_resid(sim.fm.linnorm, design_matrix_lin_norm, nsim)

for (k in 1:length(lin_norm_resid)) {
  saveRDS(lin_norm_resid[[k]], paste("Result/Lin_nonnormal/", "scenario_8_", as.character(k), ".RDS", sep = ""))
}

rm(linnorm_data_list)
rm(lin_norm_resid)



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
spec_resid <- extract_2levels_resid_spec(sim.fm.spec, design_matrix_spec, nsim)

for (k in 1:length(spec_resid)) {
  saveRDS(spec_resid[[k]], paste("Result/Special/", "scenario_9_", as.character(k), ".RDS", sep = ""))
}

rm(spec_data_list)
rm(spec_resid)


## Check setting of each saved file
## readRDS("Result/Normality/scenario_1_2.RDS")[[1]]$data_info


## How residual data are saved
## 0. Reference
## scenario_0_1.RDS -- scenario_0_6.RDS

## 1. Error_Normality
## scenario_1_1.RDS -- scenario_1_24.RDS

## 2. Error_Heteroscedasticity
## scenario_2_1.RDS -- scenario_2_18.RDS

## 3. Error_Linearity
## scenario_3_1.RDS -- scenario_3_6.RDS

## 4. RE_normality
## scenario_4_1.RDS -- scenario_4_12.RDS

## 5. Hetero & non_normal
## scenario_5_1.RDS -- scenario_5_18.RDS

## 6. Hetero & non_linear
## scenario_6_1.RDS -- scenario_6_18.RDS

## 8. Non_linear & non-normal
## scenario_8_1.RDS -- scenario_8_12.RDS

## 9. Special Cases
## scenario_9_1.RDS -- scenario_9_18.RDS



# stop cluster
stopCluster(cl)
