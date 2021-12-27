
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
library(parallel) ##< for mclapply
library(purrr)
library(readr)

# Call Functions

source("Helper_Functions.R")
source("Simulation_Functions.R")


# register cores, initialize clusters and set seed

ncores <-  12      # Use 12 core / Bio 328
cl <- makeCluster(ncores, type = "FORK")
# registerDoParallel(cl)

set.seed(4352216)
nsim <- 5      # change as needed

# Added simulations

## 1. Linearity & non_normal 
## 2. Special cases: 
#Longitudinal: 1. one variables as numeric sequence; 2. error term follows an ar(1) process within each group;
#random effect: having only random intercept)




## design_matrix_extra containing: 

# design_matrix_lin_norm
# design_matrix_spec

variance<-c("large_error","large_re") 
balance<-c("same","balanced","unbalanced")
lin_norm <- c("linear_norm", "reduced_skew", "reduced_bimodal") #mildly skewed, skewness = 1.5
special <- c("standard", "time_seq", "ar_error", "re_int")

design_matrix_spec<- expand.grid(variance, balance,special[-1], lin_norm[1])
design_matrix_lin_norm<-expand.grid(variance, balance, special[1], lin_norm[-1] )

scenario <- c("variance","balance", "special", "lin_norm")

colnames(design_matrix_spec) <- colnames (design_matrix_lin_norm) <-
  scenario

# Combine all settings considered
design_matrix_extra<- rbind(
  design_matrix_spec, 
  design_matrix_lin_norm) %>% distinct()




##lin_norm 

#get data
linnorm_data_list <- unlist(get_design_matrix_data_2(design_matrix_lin_norm), recursive = F)

sim.fm.linnorm <- mclapply (
  linnorm_data_list,
  FUN = function(x) {
    sim.hlm_linnorm(x$simdata.df, lin_norm.dsn  = as.character(x$simdata.df$lin_norm[1]), var = as.character(x$simdata.df$variance[1]))
  } , mc.cores = ncores
)


#resid

lin_norm_resid <- extract_2levels_resid(sim.fm.linnorm, design_matrix_lin_norm, nsim) 

for(k in 1:length(lin_norm_resid)){
  saveRDS(lin_norm_resid[[k]], paste("Results2/Lin_nonnormal/", "scenario_8_", as.character(k), ".RDS", sep = ""))
}         

rm(linnorm_data_list)
rm(lin_norm_resid)


## speical cases
#get data
spec_data_list <- unlist(get_design_matrix_data_2(design_matrix_spec), recursive = F)

sim.fm.spec <- mclapply(
  spec_data_list, 
  FUN = function(x) {
    sim.hlm_spec(x$simdata.df, var = as.character(x$simdata.df$variance[1]), special = as.character(x$simdata.df$special[1]), J = x$J)
  }, mc.cores = ncores
)


#get resid
spec_resid <- extract_2levels_resid_spec(sim.fm.spec, design_matrix_spec, nsim)      

for(k in 1:length(spec_resid)){
  saveRDS(spec_resid[[k]], paste("Results2/Special/", "scenario_9_", as.character(k), ".RDS", sep = ""))
}                

rm(spec_data_list)
rm(spec_resid)       



stopCluster(cl)
