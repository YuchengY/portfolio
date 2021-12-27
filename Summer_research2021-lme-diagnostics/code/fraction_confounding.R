## fraction of confounding for conditional residuals

# Generate data
set.seed(202167)
sa_simdata.df<- generate_data(50, 0, 0.5)
ba_simdata.df<- generate_data(50, 1, 0.5)
ub_simdata.df<- generate_data(50, 2, 1)

# Simulate data with that design...
sim_design <- ub_simdata.df
b0 <- 0.5
b1 <- 1
b2 <- 1
sd_error <- 1
sd_int   <- 20
sd_slope <- 20
re_df <- data.frame(group = 1:50, 
                    intercept = rnorm(50, mean = 0, sd = sd_int),
                    slope = rnorm(50, mean = 0, sd = sd_slope))

sim_df <- left_join(sim_design, re_df, by = "group") %>%
  mutate(error = rnorm(nrow(.), mean = 0, sd = sd_error),
         Y_ij = b0 + b1 * X_ij + b2 * Z_ij + intercept + slope * X_ij + error)


# Model organization
fm  <-  lmer(Y_ij ~ X_ij + Z_ij + (X_ij | group), data = sim_df)
library(mlmRev)
# fm  <- lmer(normexam ~ standLRT + sex + schgend + (1|school), Exam)
# fm <- lmer(log.radon ~ uranium + basement + house.uranium + (basement | county), radon2 )

# Extract key matrices
mod_mats  <- HLMdiag:::.lmerMod_matrices(fm)

# FC
n <- mod_mats$n
sig0 <- sigma(fm)

R <- Diagonal(n = n, x = sig0^2)
ZDZt <- sig0^2 * crossprod(lme4::getME(fm, "A"))
fc_num <- R %*% mod_mats$P %*% ZDZt %*% mod_mats$P %*% R
fc_denom <- R %*% mod_mats$P %*% R

mean(diag(fc_num) / diag(fc_denom))

