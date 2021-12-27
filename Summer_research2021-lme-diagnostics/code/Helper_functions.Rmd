# This R Script is for recording all helper functions used
# Still ongoing, update when needed


# The list of all helper function stored here:

# 1.bp
# 2.matrix_sqrt
# 3.lesaffre_verbeke
# 4.mahalanobis_ranef.lmerMod
# 5.resid_ranef.lmerMod
# 6.bound_re
# 7.plot_ranef
# 8.sim_heteroscedastic_error
# 9.bimodalDistFunc
# 10.generate_data (Update 7.26)
# 11. vec_p
# 12. vec.rpearson
# 13. bvn_intercept
# 14. bvn_slope
# 15. chi_gof
# 16. concat
# 17. generate_data_ar
# 18. generate_data_t
# 19. generate_data_intercept
# 20. vec_p2
# 21. flattenlist
# 22. get_design_matrix_data
# 23. extract_2levels_resid
# 24. all_residual_tests
# 25. test.df
# 26. rate
# 27. conduct_4_test
# 28. run_tests
# 29. get_test_results
# 30. generate_data_2 (8.3)
# 31. get_design_matrix_data_2
# 32. extract_2levels_resid_spec
# 33. conduct_4_test_mc
# 34. all_residual_tests_mc
# 35. extract_2levels_resid_omit
# 36. run_tests_new (new fun for rerun in 8.10)
# 37. get_test_results_new (new fun for rerun in 8.10)



# 1. Breusch-Pagan-Test
bp <- function(resid, fitted) {
  .sq_resid <- resid^2 # squared residual

  # Regression e^2 on the fitted values
  .mod <- lm(.sq_resid ~ fitted)

  # Regression sum of squares
  .SSR <- sum((predict(.mod) - mean(.sq_resid))^2) # anova(.mod)[1,2]
  # Error Sum of Squares
  .SSE <- sum(.sq_resid) # anova(.mod)[2,2]
  # n
  .n <- length(resid)
  # test statistics
  .x <- (.SSR / 2) / (.SSE / .n)^2

  .p <- pchisq(.x, df = 1)
  .p <- ifelse(.p >= 0.5, 1 - .p, .p)
  return(c(stats = .x, p_val = .p))
}



# 2. Matrix sqrt needed for lesaffre-verbeke stat in Singer et al.
matrix_sqrt <- function(mat) {
  mat <- as.matrix(mat)
  decomp <- eigen(mat)
  V <- decomp$vectors
  S <- sqrt(decomp$values)
  V %*% diag(S) %*% t(V)
}



# 3. Currently only for 2-level models fit via lme4
.lmerMod_matrices <- HLMdiag:::.lmerMod_matrices

lesaffre_verbeke <- function(object) {
  mats <- HLMdiag:::.lmerMod_matrices(object)
  mar_var <- mats$V - mats$X %*% tcrossprod(mats$XVXinv, mats$X)
  mar_resid <- HLMdiag:::resid_marginal(object, type = "raw")

  grp_index <- getME(object, "flist")[[1]]
  grp_names <- levels(grp_index)
  grp_sizes <- table(grp_index)

  diagnostic <- vector("numeric", length(grp_names))
  for (i in seq_along(grp_names)) {
    indx <- which(grp_index == grp_names[i])
    vi_block <- mar_var[indx, indx]

    if (grp_sizes[i] == 1) {
      ri <- solve(sqrt(vi_block)) %*% mar_resid[indx]
      Vi <- diag(grp_sizes[i]) - tcrossprod(ri)
      diagnostic[i] <- sqrt(sum(diag(tcrossprod(Vi))))
    } else {
      ri <- solve(matrix_sqrt(vi_block)) %*% mar_resid[indx]
      Vi <- diag(grp_sizes[i]) - tcrossprod(ri)
      diagnostic[i] <- sqrt(sum(diag(tcrossprod(Vi))))
    }
  }

  diagnostic / grp_sizes
}



# 4.  @importFrom diagonals split_vector fatdiag (updated at 8.3)
mahalanobis_ranef.lmerMod <- function(object) {
  ngrps <- 0
  mats <- .lmerMod_matrices(object)
  n_lev <- length(lme4::getME(object, "flist"))

  if (n_lev == 1) {
    Z <- lme4::getME(object, "Z")
    vc <- lme4::VarCorr(object)
    D <- kronecker(Diagonal(mats$ngrps), bdiag(vc))

    eblup <- tcrossprod(D, Z) %*% mats$Vinv %*% resid_marginal(object)
    # vcov_eblup <- D - tcrossprod(D, Z) %*% mats$P %*% tcrossprod(Z, D)
    vcov_eblup <- tcrossprod(D, Z) %*% mats$P %*% tcrossprod(Z, D)

    if (ncol(ranef(object)[[1]]) == 1) {
      mah_dist_eblup <- as.numeric(eblup / diag(vcov_eblup))
    } else {
      eblup_lst <- diagonals::split_vector(eblup, size = 2)
      vcov_eblup_lst <- diagonals::fatdiag(vcov_eblup, steps = ngrps(object)) %>%
        diagonals::split_vector(size = 4) %>%
        map(~ matrix(.x, nrow = 2, byrow = TRUE))
      mah_dist_eblup <- purrr::map2_dbl(eblup_lst, vcov_eblup_lst, ~ t(.x) %*% MASS::ginv(.y) %*% .x)
    }
  }
  return(mah_dist_eblup)
}



# 5. extract random effects
resid_ranef.lmerMod <- function(object, whichel, standardize = FALSE) {
  if (!is.logical(standardize)) {
    stop("standardize must be logical (TRUE or FALSE).")
  }

  flist <- lme4::getME(object, "flist")

  if (standardize) {
    re <- lme4::ranef(object, condVar = TRUE)
    vc <- lme4::VarCorr(object)
    for (i in names(flist)) {
      diag_var <- diag(as.matrix(vc[[i]])) - diag(as.matrix(attr(re[[i]], "postVar")[, , 1]))
      re[[i]] <- sweep(re[[i]], 2, sqrt(diag_var), FUN = "/")
    }
  } else {
    re <- lme4::ranef(object)
  }

  re[[whichel]]
}



# 6. Simultaneous bands of Schuetzenmesiter and Piepho,
# Use STB::getSTB().... this is the general idea...
bound_re <- function(.merMod, whichel, standardize = FALSE, B = 2000) {
  orig_re <- resid_ranef.lmerMod(.merMod, whichel = whichel, standardize = standardize)

  boot_re <- lmeresampler::parametric_bootstrap(
    .merMod,
    .f = function(x) resid_ranef.lmerMod(x, whichel = whichel, standardize = standardize),
    B = B
  )
  n_re <- ncol(boot_re$replicates) - 1
  res <- mutate(boot_re$replicates,
    row_names = str_extract(rownames(boot_re$replicates), "\\d+")
  )

  sim_mats <- vector(mode = "list", length = n_re)
  for (i in seq_len(n_re)) {
    sim_mats[[i]] <- res %>%
      dplyr::select(colnames(res)[i], row_names) %>%
      group_by(row_names) %>%
      mutate(row = row_number()) %>%
      pivot_wider(names_from = "row_names", values_from = colnames(res)[i]) %>%
      dplyr::select(-row)
  }

  bounds <- lapply(sim_mats,
    FUN = function(.x) {
      mat <- t(apply(.x, 1, sort))
      STB::getSTB(mat, output = FALSE)
    }
  )

  result <- lapply(bounds, FUN = function(.y) as_tibble(t(.y$Q)) %>% set_names(nm = c("lower", "upper")))
  result <- map2(result, orig_re, .f = ~ bind_cols(.x, pred = sort(.y)))
  names(result) <- colnames(res)[seq_len(n_re)]

  bind_rows(result, .id = "re")
}



# 7. Plots of the random effects
plot_ranef <- function(.merMod, whichel = 1, standardize = FALSE, B = 2000, ncol = 2) {
  re_df <- bound_re(.merMod, whichel, standardize, B) %>%
    group_by(re) %>%
    mutate(expected = qqnorm(pred, plot.it = F)$x)

  ggplot(data = re_df, aes(sample = pred)) +
    qqplotr::stat_qq_line(distribution = "norm", dparams = list(mean = 0, sd = 1)) +
    geom_ribbon(aes(x = expected, y = pred, ymin = lower, ymax = upper), alpha = 0.3) +
    # stat_qq_band(bandType = "ts", fill = "orange", alpha = 0.3, distribution = "norm", dparams = list(mean = 0, sd = 1)) +
    qqplotr::stat_qq_point(distribution = "norm", dparams = list(mean = 0, sd = 1)) +
    labs(x = "N(0, 1) quantiles", y = "Sample data") +
    facet_wrap(~re, ncol = ncol)
}



# 8. simulate heteroscedastic errors
sim_heteroscedastic_error <- function(x, var, hfactor) {
  var_heter <- var + hfactor * (x - min(x))
  rnorm(length(x), mean = 0, sd = sqrt(var_heter))
}



# 9. generate bimodal distribution
bimodalDistFunc <- function(n, cpct, mu1, mu2, sig1, sig2) {
  y0 <- rnorm(n, mean = mu1, sd = sig1)
  y1 <- rnorm(n, mean = mu2, sd = sig2)
  flag <- rbinom(n, size = 1, prob = cpct)
  y <- y0 * (1 - flag) + y1 * flag
}



# 10. generate two-level data set (Update 7.26)
generate_data <- function(ngroup = 50, variance = "large_error", balance = "same",
                          norm = "norm", hetero = 0, lin = "full", omit = "full", RE = "norm_re", hetero_lin = "linear_homo", hetero_norm = "homo_normal") {
  N <- ngroup # default 50 groups
  # Create two random effects (slope and intercept)
  sigma_1 <- 0.5
  sigma_2 <- 0.5
  rho <- 0.5
  cov.matrix <- matrix(c(sigma_1^2, sigma_1 * sigma_2 * rho, sigma_1 * sigma_2 * rho, sigma_2^2),
    nrow = 2,
    byrow = TRUE
  )
  random.effects <- mvrnorm(N, mu = c(0, 0), Sigma = cov.matrix)

  # Set up Groups
  cluster.df <- data.frame(group = c(1:N), a = rnorm(N)) # Each group has its "characteristic" (unused now)

  cluster.df <- within(cluster.df, {
    slope_x <- 1
    slope_z <- 1
    intercept <- 5
  })

  cluster.df$slope_x <- cluster.df$slope_x + random.effects[, 1]
  cluster.df$intercept <- cluster.df$intercept + random.effects[, 2]

  # Set up number of observations within each groups

  if (balance == "same") {
    J <- rep(25, ngroup)
  }

  if (balance == "balanced") {
    J <- sample(20:30, ngroup, replace = TRUE)
  }

  if (balance == "unbalanced") {
    J <- sample(2:50, ngroup, replace = TRUE)
  }

  M <- sum(J) # Total number of observations

  x.grid <- list()
  z.grid <- list()
  j <- list()
  group <- list()

  for (i in 1:N) {
    x.grid[[i]] <- (rnorm(J[[i]]))
    z.grid[[i]] <- abs(rnorm(J[[i]]))
    j[[i]] <- rep(1:J[[i]])
    group[[i]] <- rep(c(i), J[[i]])
  }

  within.cluster.df <- data.frame(group = unlist(group), j = unlist(j), X_ij = unlist(x.grid), Z_ij = unlist(z.grid))

  # Without Y_ij
  simdata.df <- merge(cluster.df, within.cluster.df)
  # norm = "norm"; hetero = 0; lin = "full"; omit ="full"; RE = "norm_re"
  # hetero = 0, lin = "full", omit ="full", RE = "norm_re"
  simdata.df <- simdata.df %>% mutate(variance = variance, balance = balance, normality = norm, hetero = hetero, linearity = lin, omit = omit, random_effects = RE, hetero_lin = hetero_lin, hetero_norm = hetero_norm)

  simdata.df <- as_tibble(simdata.df[, c("group", "j", "X_ij", "Z_ij", "variance", "balance", "normality", "hetero", "linearity", "omit", "random_effects", "hetero_lin", "hetero_norm")])

  return(simdata.df)
}


# 11. vectorize all moments according to the random generated mean
vec_p <- function(n, df) {
  m_list <- vector("list", length = n)
  for (i in 1:n) {
    m_list[[i]] <- c(df[i], variance = df[n + 1], skewness = df[n + 2], kurtosis = df[n + 3])
  }
  m_list
}


# 12. vectorize rpearson
vec.rpearson <- Vectorize(rpearson)



# 13. generate bivariate normal with skewed intercept
bvn_intercept <- function(n, m1, s1, m2, s2, rho) {
  X2 <- rnorm(n, m2, s2)
  p_df <- c(mean = m1 + (s2 / s1) * rho * (X2 - m1), variance = (1 - rho^2) * s1^2, skewness = 1.5, kurtosis = 6)
  p <- vec_p(n, p_df)
  X1 <- vec.rpearson(1, moments = p)

  cbind(X1, X2)
}



# 14. generate bivariate normal with skewed slope

bvn_slope <- function(n, m1, s1, m2, s2, rho) {
  X1 <- rnorm(n, m1, s1)

  p_df <- c(mean = m2 + (s2 / s1) * rho * (X1 - m1), variance = (1 - rho^2) * s2^2, skewness = 1.5, kurtosis = 6)
  p <- vec_p(n, p_df)
  X2 <- vec.rpearson(1, moments = p)

  cbind(X1, X2)
}



# 15. setup function for gof test
chi_gof <- function(x) {
  a <- mahalanobis_ranef.lmerMod(x)
  a
}



# 16. concatenate helper function #16 generate_data_ar
concat <- function(n, df) {
  e_df <- as.matrix(df[, 1], ncol = 1)
  for (i in 2:n) {
    e_df <- rbind(e_df, as.matrix(df[, i], ncol = 1))
  }
  e_df
}



# 17. generating AR(1) error two-level model (modified from # 10 generate_data function from Yicheng, can merge into one if needed)
generate_data_ar <- function(ngroup = 50, balance = 0, rho = 0.5) {
  N <- ngroup # default 50 groups

  # Create two random effects (slope and intercept)
  sigma_1 <- 0.2
  sigma_2 <- 0.5
  rho <- rho
  cov.matrix <- matrix(c(sigma_1^2, sigma_1 * sigma_2 * rho, sigma_1 * sigma_2 * rho, sigma_2^2),
    nrow = 2,
    byrow = TRUE
  )
  random.effects <- mvrnorm(N, mu = c(0, 0), Sigma = cov.matrix)

  # Set up Groups
  cluster.df <- data.frame(group = c(1:N), a = rnorm(N)) # Each group has its "characteristic" from normal distribution
  # cluster.df

  cluster.df <- within(cluster.df, {
    slope_x <- 0.5 * a + 1
    intercept <- 2 * a + 5
  })

  cluster.df$slope_x <- cluster.df$slope_x + random.effects[, 1]
  cluster.df$intercept <- cluster.df$intercept + random.effects[, 2]

  # Set up number of observations within each groups

  if (balance == 0) {
    J <- rep(25, ngroup)
  }

  if (balance == 1) {
    J <- sample(20:30, ngroup, replace = TRUE)
  }

  if (balance == 2) {
    J <- sample(2:50, ngroup, replace = TRUE)
  }

  M <- sum(J) # Total number of observations

  x.grid <- list()
  z.grid <- list()
  j <- list()
  group <- list()

  for (i in 1:ngroup) {
    x.grid[[i]] <- (rnorm(J[[i]]))
    z.grid[[i]] <- abs(rnorm(J[[i]]))
    j[[i]] <- rep(1:J[[i]])
    group[[i]] <- rep(c(i), J[[i]])
  }

  slope_z <- 0.5

  within.cluster.df <- data.frame(group = unlist(group), j = unlist(j), X_ij = unlist(x.grid), Z_ij = unlist(z.grid))

  simdata.df <- merge(cluster.df, within.cluster.df)

  # ar correlation
  ar.val <- 0.4
  # error sd
  sigma <- 1.5
  # ar error by group
  ### note: use arima.sim(model=list(ar=ar.val), n=x) * sqrt(1-ar.val^2) * sigma
  ### construction, so that the true error SD is equal to sigma
  eij <- unlist(sapply(J, function(x) arima.sim(model = list(order = c(1, 0, 0), ar = ar.val), n = x) * sqrt(1 - ar.val^2) * sigma))
  # concat eij
  e <- concat(N, eij)

  simdata.df <- within(simdata.df, Y_ij <- intercept + X_ij * slope_x + Z_ij * slope_z + e)

  simdata.df <- simdata.df[, c("group", "j", "X_ij", "Z_ij", "Y_ij")]

  return(as_tibble(simdata.df))
}



# 18. generating longitidinal two-level model (one more fixed effect as a numerical sequence for each group)
generate_data_t <- function(ngroup = 50, balance = 0, rho12 = 0.5, rho13 = 0.4, rho23 = 0.6) {
  N <- ngroup # default 50 groups

  # Create two random effects (slope and intercept)
  sigma_1 <- 0.2
  sigma_2 <- 0.6
  sigma_3 <- 0.3
  rho12 <- rho12
  rho13 <- rho13
  rho23 <- rho23

  # three fixed effect terms, have 3*3 covariance matrix
  cov.matrix <- matrix(c(sigma_1^2, sigma_1 * sigma_2 * rho12, sigma_1 * sigma_3 * rho13, sigma_1 * sigma_2 * rho12, sigma_2^2, sigma_2 * sigma_3 * rho23, sigma_1 * sigma_3 * rho13, sigma_2 * sigma_3 * rho23, sigma_3^2), nrow = 3, byrow = TRUE)
  random.effects <- mvrnorm(N, mu = c(0, 0, 0), Sigma = cov.matrix)

  # Set up Groups, add in time variable
  cluster.df <- data.frame(group = c(1:N), a = rnorm(N), t = 1:N) # Each group has its "characteristic" from normal distribution
  # cluster.df

  cluster.df <- within(cluster.df, {
    slope_x <- 0.5 * a + 1
    # slope for time variable
    slope_t <- 0.5 * t + 1
    intercept <- 2 * a + 5
  })

  cluster.df$slope_x <- cluster.df$slope_x + random.effects[, 1]
  cluster.df$slope_t <- cluster.df$slope_t + random.effects[, 3]
  cluster.df$intercept <- cluster.df$intercept + random.effects[, 2]

  # Set up number of observations within each groups

  if (balance == 0) {
    J <- rep(25, ngroup)
  }

  if (balance == 1) {
    J <- sample(20:30, ngroup, replace = TRUE)
  }

  if (balance == 2) {
    J <- sample(2:50, ngroup, replace = TRUE)
  }

  M <- sum(J) # Total number of observations

  x.grid <- list()
  z.grid <- list()
  t.grid <- list()
  j <- list()
  group <- list()

  for (i in 1:ngroup) {
    x.grid[[i]] <- (rnorm(J[[i]]))
    t.grid[[i]] <- (1:N)
    z.grid[[i]] <- abs(rnorm(J[[i]]))
    j[[i]] <- rep(1:J[[i]])
    group[[i]] <- rep(c(i), J[[i]])
  }

  slope_z <- 0.5
  # slope for time variable
  slope_t <- 0.5

  within.cluster.df <- data.frame(group = unlist(group), j = unlist(j), X_ij = unlist(x.grid), T_ij = unlist(t.grid), Z_ij = unlist(z.grid))

  simdata.df <- merge(cluster.df, within.cluster.df)

  simdata.df <- within(simdata.df, Y_ij <- intercept + X_ij * slope_x + T_ij * slope_t + Z_ij * slope_z + rnorm(n = M))

  simdata.df <- simdata.df[, c("group", "j", "X_ij", "T_ij", "Z_ij", "Y_ij")]

  return(as_tibble(simdata.df))
}



# 19. generating two-level dataset with only random intercept
generate_data_intercept <- function(ngroup = 50, balance = 0) {
  N <- ngroup # default 50 groups

  # Create two random effects (slope and intercept)
  sigma_1 <- 0.2 # sd of the random effect (intercept)
  random.effects <- rnorm(N, mean = 0, sd = sigma_1)

  # Set up Groups
  cluster.df <- data.frame(group = c(1:N), a = rnorm(N)) # Each group has its "characteristic" from normal distribution
  # cluster.df

  cluster.df <- within(cluster.df, {
    slope_x <- 0.5 * a + 1
    intercept <- 2 * a + 5
  })

  # cluster.df$slope_x = cluster.df$slope_x + random.effects[, 1]
  cluster.df$intercept <- cluster.df$intercept + random.effects

  # Set up number of observations within each groups

  if (balance == 0) {
    J <- rep(25, ngroup)
  }

  if (balance == 1) {
    J <- sample(20:30, ngroup, replace = TRUE)
  }

  if (balance == 2) {
    J <- sample(2:50, ngroup, replace = TRUE)
  }

  M <- sum(J) # Total number of observations

  x.grid <- list()
  z.grid <- list()
  j <- list()
  group <- list()

  for (i in 1:ngroup) {
    x.grid[[i]] <- (rnorm(J[[i]]))
    z.grid[[i]] <- abs(rnorm(J[[i]]))
    j[[i]] <- rep(1:J[[i]])
    group[[i]] <- rep(c(i), J[[i]])
  }

  # slope_z <- 0.5

  within.cluster.df <- data.frame(group = unlist(group), j = unlist(j), X_ij = unlist(x.grid), Z_ij = unlist(z.grid))

  simdata.df <- merge(cluster.df, within.cluster.df)

  simdata.df <- within(simdata.df, Y_ij <- intercept + X_ij * slope_x + Z_ij + rnorm(n = M))

  simdata.df <- simdata.df[, c("group", "j", "X_ij", "Z_ij", "Y_ij")]

  return(as_tibble(simdata.df))
}



# 20. vectorize all moments for different covariance
vec_p2 <- function(n, df) {
  m_list <- vector("list", length = n)
  for (i in 1:n) {
    m_list[[i]] <- c(df[1], variance = df[i + 1], skewness = 1.5, kurtosis = 6)
  }
  m_list
}



# 21. flatten nested list
# source: https://stackoverflow.com/questions/16300344/how-to-flatten-a-list-of-lists
flattenlist <- function(x) {
  morelists <- sapply(x, function(xprime) class(xprime)[1] == "list")
  out <- c(x[!morelists], unlist(x[morelists], recursive = FALSE))
  if (sum(morelists)) {
    Recall(out)
  } else {
    return(out)
  }
}



# 22. get data based on specific sub design matrix
get_design_matrix_data <- function(design_matrix_input) {
  data_list <- vector("list", length = nrow(design_matrix_input))
  for (i in 1:nrow(design_matrix_input)) {
    data_list[[i]] <- vector("list", length = nsim)
    data_list[[i]] <- replicate(
      n = nsim,
      generate_data(
        ngroup = 50,
        variance = design_matrix_input %>% pull(variance) %>% `[`(i),
        balance = design_matrix_input %>% pull(balance) %>% `[`(i),
        norm = design_matrix_input %>% pull(normality) %>% `[`(i),
        hetero = design_matrix_input %>% pull(heteroscedasticity) %>% `[`(i),
        lin = design_matrix_input %>% pull(linearity) %>% `[`(i),
        omit = design_matrix_input %>% pull(fixed_effect) %>% `[`(i),
        RE = design_matrix_input %>% pull(random_effect) %>% `[`(i),
        hetero_lin = design_matrix_input %>% pull(Hetero_lin) %>% `[`(i),
        hetero_norm = design_matrix_input %>% pull(Hetero_norm) %>% `[`(i)
      ),
      simplify = FALSE
    )
  }
  return(data_list)
}



# 23. A function for extracting residuals from simulated models
extract_2levels_resid <- function(sim.fm, design_matrix, nsim) {
  # note: here we need to pass in the sub design matrix for that specific scenario, eg. design_matrix_hetero_norm

  model_resid <- mclapply(
    sim.fm,
    # first loop through all design settings of that scenario
    # foreach(i = 1: (length(sim.fm)/nsim)) %:%
    # then loop through the data generated under that design setting
    FUN = function(x) {
      # data_info<-design_matrix[i,]
      sim_model <- lmer(Y_ij ~ X_ij + Z_ij + (X_ij | group),
        data = as_tibble(x),
        control = lmerControl(calc.derivs = FALSE)
      )

      # level 1 resid
      first_resid <- merge(
        hlm_resid(sim_model, standardize = T, include.ls = T),
        hlm_resid(sim_model, include.ls = T)
      )
      # level 2 resid
      sec_resid <- merge(
        hlm_resid(sim_model, level = "group", standardize = T, include.ls = T),
        hlm_resid(sim_model, level = "group", include.ls = T)
      )

      # Add mahalanobis distance to level 2
      sec_resid$mahalanobis <- mahalanobis_ranef.lmerMod(sim_model)

      # combine into one model_resid object
      model_resid <- list(first_resid, sec_resid) # list(data_info, first_resid, sec_resid)
      names(model_resid) <- c("first_resid", "second_resid") # c("data_info", "first_resid", "second_resid")
      model_resid
    },
    mc.cores = ncores
  )

  # Splitting into list by design
  split_resid <- split(model_resid, rep(seq(nrow(design_matrix)), each = nsim))
  res <- lapply(
    1:nrow(design_matrix),
    function(i) {
      lapply(split_resid[[i]], function(.x) c("data_info" = list(design_matrix[i, ]), .x))
    }
  )

  return(res)
}



# 24. conduct shapiro, bp, anova, and normality test of re
all_residual_tests <- function(model_residual, design_matrix) {
  # norm
  tests <-
    foreach(j = 1:nrow(design_matrix)) %:%
    foreach(i = 1:nsim, .combine = c) %dopar% { # add .combine = c
      setting <- model_residual[[j]]

      # norm
      shapiro <- shapiro.test(setting[[i]]$first_resid$.std.resid)$p.value

      # hetero
      bp_test <- bp(setting[[i]]$first_resid$.resid, setting[[i]]$first_resid$.fitted)["p_val"]

      # re_norm
      re_norm <- ks.test(setting[[i]]$second_resid$mahalanobis, "pchisq", df = 2)$p.value

      # linearity
      chol.mar.resid <- setting[[i]]$first_resid$.chol.mar.resid
      fittted <- setting[[i]]$first_resid$.mar.fitted
      lin_test <- anova(
        lm(chol.mar.resid ~ fittted),
        lm(chol.mar.resid ~ fittted + I(fittted^2))
      )$`Pr(>F)`[2]

      # bind them together and give names
      all_tests <- list(shapiro, bp_test, re_norm, lin_test)
      names(all_tests) <- c("shapiro", "bp_test", "re_norm", "lin_test")
      all_tests
    }

  nms <- c("shapiro", "bp_test.p_val", "re_norm", "lin_test")

  reject_rate <-
    foreach(j = 1:nrow(design_matrix)) %:%
    foreach(i = 1:length(nms), .combine = c) %dopar% {
      mean(unlist(tests[[j]])[names(unlist(tests[[j]])) == as.character(nms[i])] < 0.05)
    }
  names(reject_rate) <- c("shapiro", "bp_test", "re_norm", "lin_test")

  return(reject_rate)
}



# 25. test results df ready for adding into design matrix
test.df <- function(test_results, design_matrix) {
  nms <- c("shapiro", "bp_test.p_val", "re_norm", "lin_test")
  df <- foreach(i = 1:length(nms)) %:%
    foreach(j = 1:nrow(design_matrix)) %dopar% {
      test_results[[j]][i]
    }

  return(df)
}



# 26. rate of p value < 0.05
rate <- function(df) {
  mean(df < 0.05)
}



# 27. # conduct four test and save with design matrix
conduct_4_test <- function(model_residuals, design_matrix) {
  test_results <- all_residual_tests(model_residuals, design_matrix)

  df <- test.df(test_results, design_matrix)

  design_results <- design_matrix
  design_results$shapiro <- as.numeric(df[[1]])
  design_results$bp_test <- as.numeric(df[[2]])
  design_results$re_norm <- as.numeric(df[[3]])
  design_results$lin_test <- as.numeric(df[[4]])

  return(design_results)
}


# 28. run test using non-nested mclapply
run_tests <- function(model_residual, ncores) {
  setting <- unique(lapply(model_residual, `[[`, 1))[[1]]

  tests <- mclapply(
    model_residual,
    function(x) {
      # norm
      shapiro <- shapiro.test(x$first_resid$.std.resid)$p.value

      # hetero
      bp_test <- bp(x$first_resid$.resid, x$first_resid$.fitted)["p_val"]

      # re_norm
      if (setting$special != "re_int") {
        re_norm <- ks.test(x$second_resid$mahalanobis, "pchisq", df = 2)$p.value
      }
      if (setting$special == "re_int") {
        re_norm <- ks.test(x$second_resid$mahalanobis, "pchisq", df = 1)$p.value
      }

      # linearity
      chol.mar.resid <- x$first_resid$.chol.mar.resid
      fittted <- x$first_resid$.mar.fitted
      lin_test <- anova(
        lm(chol.mar.resid ~ fittted),
        lm(chol.mar.resid ~ fittted + I(fittted^2))
      )$`Pr(>F)`[2]

      # bind them together and give names
      all_tests <- list(shapiro, bp_test, re_norm, lin_test)
      names(all_tests) <- c("shapiro", "bp_test", "re_norm", "lin_test")
      all_tests
    },
    mc.cores = ncores
  )

  test_df <- do.call("rbind", tests)
  results <- apply(test_df, 2, function(x) mean(x < 0.05, na.rm = TRUE))

  cbind(setting, t(results))
}



# 29. execute the analysis using run_tests
get_test_results <- function(folder_path, ncores) {
  files <- dir(folder_path)
  map_dfr(files, ~ read_rds(paste0(folder_path, "/", .x)) %>% run_tests(ncores))
}



# 30. generate data function modified for "extra" scenarios
generate_data_2 <- function(ngroup = 50, variance = "large_error", balance = "same", special = "standard", lin_norm = "lin_norm") {
  N <- ngroup # default 50 groups
  # Create two random effects (slope and intercept)
  sigma_1 <- 0.5
  sigma_2 <- 0.5
  rho <- 0.5

  if (special != "re_int") {
    cov.matrix <- matrix(c(sigma_1^2, sigma_1 * sigma_2 * rho, sigma_1 * sigma_2 * rho, sigma_2^2),
      nrow = 2,
      byrow = TRUE
    )
    random.effects <- mvrnorm(N, mu = c(0, 0), Sigma = cov.matrix)

    # Set up Groups
    cluster.df <- data.frame(group = c(1:N), a = rnorm(N)) # Each group has its "characteristic" (unused now)

    cluster.df <- within(cluster.df, {
      slope_x <- 1
      intercept <- 5
    })

    cluster.df$slope_x <- cluster.df$slope_x + random.effects[, 1]
    cluster.df$intercept <- cluster.df$intercept + random.effects[, 2]
  }
  if (special == "re_int") {
    random.effects <- rnorm(N, 0, sd = 0.5)

    cluster.df <- data.frame(group = c(1:N), a = rnorm(N))

    cluster.df <- within(cluster.df, {
      slope_x <- 1
      slope_z <- 1
      intercept <- 5
    })

    cluster.df$slope_x <- cluster.df$slope_x
    cluster.df$intercept <- cluster.df$intercept + random.effects
  }

  if (special != "re_int") {

    # Set up number of observations within each groups

    # large datasets, no dorpout
    if (balance == "same") {
      J <- rep(25, ngroup)
    }
    # small datasets, balanced
    if (balance == "balanced") {
      N <- 20
      J <- sample(8:10, N, replace = TRUE)
    }
    # with dropouts
    if (balance == "unbalanced") {
      N <- 20
      J <- sample(1:10, N, replace = TRUE)
    }
  } else {
    if (balance == "same") {
      J <- rep(25, ngroup)
    }

    if (balance == "balanced") {
      J <- sample(20:30, ngroup, replace = TRUE)
    }

    if (balance == "unbalanced") {
      J <- sample(2:50, ngroup, replace = TRUE)
    }
  }

  M <- sum(J) # Total number of observations

  x.grid <- list()
  z.grid <- list()
  j <- list()
  group <- list()

  if (special == "time_seq") {
    t.grid <- list()

    for (i in 1:N) {
      x.grid[[i]] <- (rnorm(J[[i]]))
      t.grid[[i]] <- rep(1:J[[i]])
      z.grid[[i]] <- abs(rnorm(J[[i]]))
      j[[i]] <- rep(1:J[[i]])
      group[[i]] <- rep(c(i), J[[i]])
    }
  }

  # if (special == "standard"| special == "ar_error"){
  if (special != "time_seq") {
    for (i in 1:N) {
      x.grid[[i]] <- (rnorm(J[[i]]))
      z.grid[[i]] <- abs(rnorm(J[[i]]))
      j[[i]] <- rep(1:J[[i]])
      group[[i]] <- rep(c(i), J[[i]])
    }

    within.cluster.df <- data.frame(group = unlist(group), j = unlist(j), X_ij = unlist(x.grid), Z_ij = unlist(z.grid))
  }

  if (special == "time_seq") {
    x.grid <- list()
    z.grid <- list()
    t.grid <- list()
    j <- list()
    group <- list()

    for (i in 1:N) {
      x.grid[[i]] <- (rnorm(J[[i]]))
      t.grid[[i]] <- rep(1:J[[i]])
      z.grid[[i]] <- abs(rnorm(J[[i]]))
      j[[i]] <- rep(1:J[[i]])
      group[[i]] <- rep(c(i), J[[i]])
    }

    slope_z <- 0.5
    # slope for time variable
    slope_t <- 0.5

    within.cluster.df <- data.frame(group = unlist(group), j = unlist(j), X_ij = unlist(x.grid), T_ij = unlist(t.grid), Z_ij = unlist(z.grid))
  }
  # Without Y_ij
  simdata.df <- merge(cluster.df, within.cluster.df)


  simdata.df <- simdata.df %>% mutate(variance = variance, balance = balance, special = special, lin_norm = lin_norm)

  if (special == "time_seq") {
    simdata.df <- as_tibble(simdata.df[, c("group", "j", "X_ij", "T_ij", "Z_ij", "variance", "balance", "special", "lin_norm")])
  }

  if (special != "time_seq") {
    simdata.df <- as_tibble(simdata.df[, c("group", "j", "X_ij", "Z_ij", "variance", "balance", "special", "lin_norm")])
  }

  op <- list(simdata.df, J)
  names(op) <- c("simdata.df", "J")

  return(op)
}



# 31. updated for "extra" scenarios
get_design_matrix_data_2 <- function(design_matrix_input) {
  data_list <- vector("list", length = nrow(design_matrix_input))
  for (i in 1:nrow(design_matrix_input)) {
    data_list[[i]] <- vector("list", length = nsim)

    data_list[[i]] <- replicate(
      n = nsim,
      generate_data_2(
        ngroup = 50,
        variance = design_matrix_input %>% pull(variance) %>% `[`(i),
        balance = design_matrix_input %>% pull(balance) %>% `[`(i),
        special = design_matrix_input %>% pull(special) %>% `[`(i),
        lin_norm = design_matrix_input %>% pull(lin_norm) %>% `[`(i)
      ),
      simplify = FALSE
    )
  }
  return(data_list)
}



# 32. updated for "extra" scenarios
extract_2levels_resid_spec <- function(sim.fm, design_matrix, nsim) {
  # note: here we need to pass in the sub design matrix for that specific scenario, eg. design_matrix_hetero_norm

  model_resid <- lapply(
    sim.fm,
    # first loop through all design settings of that scenario

    # then loop through the data generated under that design setting
    FUN = function(x) {
      if (ncol(as_tibble(x)) == 12) {
        # data_info<-design_matrix[i,]
        sim_model <- lmer(Y_ij ~ X_ij + Z_ij + T_ij + (X_ij | group),
          data = as_tibble(x),
          control = lmerControl(calc.derivs = FALSE)
        )
      }
      if (ncol(as_tibble(x)) == 11) {
        # data_info<-design_matrix[i,]
        sim_model <- lmer(Y_ij ~ X_ij + Z_ij + (X_ij | group),
          data = as_tibble(x),
          control = lmerControl(calc.derivs = FALSE)
        )
      }
      if (ncol(as_tibble(x)) == 10) {
        # data_info<-design_matrix[i,]
        sim_model <- lmer(Y_ij ~ X_ij + Z_ij + (1 | group),
          data = as_tibble(x),
          control = lmerControl(calc.derivs = FALSE)
        )
      }
      # level 1 resid
      first_resid <- merge(
        hlm_resid(sim_model, standardize = T, include.ls = T),
        hlm_resid(sim_model, include.ls = T)
      )
      # level 2 resid
      sec_resid <- merge(
        hlm_resid(sim_model, level = "group", standardize = T, include.ls = T),
        hlm_resid(sim_model, level = "group", include.ls = T)
      )


      # Add mahalanobis distance to level 2
      sec_resid$mahalanobis <- mahalanobis_ranef.lmerMod(sim_model)

      # combine into one model_resid object
      model_resid <- list(first_resid, sec_resid) # list(data_info, first_resid, sec_resid)
      names(model_resid) <- c("first_resid", "second_resid") # c("data_info", "first_resid", "second_resid")
      model_resid
    }
  )

  # Splitting into list by design

  split_resid <- split(model_resid, rep(seq(nrow(design_matrix)), each = nsim))

  res <- lapply(
    1:nrow(design_matrix),
    function(i) {
      lapply(split_resid[[i]], function(.x) c("data_info" = list(design_matrix[i, ]), .x))
    }
  )

  return(res)
}



# 33. conduct_4_test_mc
conduct_4_test_mc <- function(model_residuals, design_matrix) {
  test_results <- all_residual_tests_mc(model_residuals, design_matrix)

  df <- test.df(test_results, design_matrix)

  design_results <- design_matrix
  design_results$shapiro <- as.numeric(flatten(df[[1]]))
  design_results$bp_test <- as.numeric(flatten(df[[2]]))
  design_results$re_norm <- as.numeric(flatten(df[[3]]))
  design_results$lin_test <- as.numeric(flatten(df[[4]]))

  return(design_results)
}



# 34. all_residual_tests_mc
all_residual_tests_mc <- function(model_residual, design_matrix) {
  # norm
  tests <-
    mclapply(1:nrow(design_matrix), function(j) {
      mclapply(1:nsim, function(i) {
        setting <- model_residual[[j]]

        # norm
        shapiro <- shapiro.test(setting[[i]]$first_resid$.std.resid)$p.value

        # hetero
        bp_test <- bp(setting[[i]]$first_resid$.resid, setting[[i]]$first_resid$.fitted)["p_val"]

        # re_norm
        re_norm <- ks.test(setting[[i]]$second_resid$mahalanobis, "pchisq", df = 2)$p.value

        # linearity
        chol.mar.resid <- setting[[i]]$first_resid$.chol.mar.resid
        fittted <- setting[[i]]$first_resid$.mar.fitted
        lin_test <- anova(
          lm(chol.mar.resid ~ fittted),
          lm(chol.mar.resid ~ fittted + I(fittted^2))
        )$`Pr(>F)`[2]

        # bind them together and give names
        all_tests <- list(shapiro, bp_test, re_norm, lin_test)
        names(all_tests) <- c("shapiro", "bp_test", "re_norm", "lin_test")
        all_tests
      }, mc.cores = ncores)
    }, mc.cores = ncores)

  nms <- c("shapiro", "bp_test.p_val", "re_norm", "lin_test")

  reject_rate <-
    mclapply(1:nrow(design_matrix), function(j) {
      mclapply(1:length(nms), function(i) {
        mean(unlist(tests[[j]])[names(unlist(tests[[j]])) == as.character(nms[i])] < 0.05)
      }, mc.cores = ncores)
    }, mc.cores = ncores)

  return(reject_rate)
}


# 35. extract_2levels_resid_omit
extract_2levels_resid_omit <- function(sim.fm, design_matrix, nsim) {
  # note: here we need to pass in the sub design matrix for that specific scenario, eg. design_matrix_hetero_norm

  model_resid <- mclapply(
    sim.fm,
    FUN = function(x) {
      # data_info<-design_matrix[i,]
      sim_model <- lmer(Y_ij ~ X_ij + (X_ij | group),
        data = as_tibble(x),
        control = lmerControl(calc.derivs = FALSE)
      )

      # level 1 resid
      first_resid <- merge(
        hlm_resid(sim_model, standardize = T, include.ls = T),
        hlm_resid(sim_model, include.ls = T)
      )
      # level 2 resid
      sec_resid <- merge(
        hlm_resid(sim_model, level = "group", standardize = T, include.ls = T),
        hlm_resid(sim_model, level = "group", include.ls = T)
      )

      # Add mahalanobis distance to level 2
      sec_resid$mahalanobis <- mahalanobis_ranef.lmerMod(sim_model)

      # combine into one model_resid object
      model_resid <- list(first_resid, sec_resid) # list(data_info, first_resid, sec_resid)
      names(model_resid) <- c("first_resid", "second_resid") # c("data_info", "first_resid", "second_resid")
      model_resid
    },
    mc.cores = ncores
  )

  # Splitting into list by design
  split_resid <- split(model_resid, rep(seq(nrow(design_matrix)), each = nsim))
  res <- lapply(
    1:nrow(design_matrix),
    function(i) {
      lapply(split_resid[[i]], function(.x) c("data_info" = list(design_matrix[i, ]), .x))
    }
  )

  return(res)
}



### 36. new fun for rerun in 8.10
run_tests_new <- function(model_residual, ncores) {
  setting <- unique(lapply(model_residual, `[[`, 1))[[1]]

  tests <- mclapply(
    model_residual,
    function(x) {
      # norm
      shapiro <- shapiro.test(x$first_resid$.std.resid)$p.value

      # hetero
      bp_test <- bp(x$first_resid$.resid, x$first_resid$.fitted)["p_val"]


      # ls bp test
      # bp_test.ls<-bp(x$first_resid$.ls.resid, x$first_resid$.ls.fitted)["p_val"]

      if (setting$special != "re_int") {
        # re_norm
        re_norm <- ks.test(x$second_resid$mahalanobis, "pchisq", df = 2)$p.value
        re_norm_int <- NA
      }
      if (setting$special == "re_int") {
        re_norm <- ks.test(x$second_resid$mahalanobis, "pchisq", df = 1)$p.value
        re_norm_int <- shapiro.test(x$second_resid$mahalanobis)$p.value
      }
      # linearity
      chol.mar.resid <- x$first_resid$.chol.mar.resid
      fittted <- x$first_resid$.mar.fitted
      lin_test <- anova(
        lm(chol.mar.resid ~ fittted),
        lm(chol.mar.resid ~ fittted + I(fittted^2))
      )$`Pr(>F)`[2]

      # bind them together and give names
      all_tests <- list(shapiro, bp_test, re_norm, re_norm_int, lin_test)
      names(all_tests) <- c("shapiro", "bp_test", "re_norm", "re_norm_int", "lin_test")
      all_tests
    },
    mc.cores = ncores
  )

  test_df <- do.call("rbind", tests)
  results <- apply(test_df, 2, function(x) mean(x < 0.05, na.rm = TRUE))

  cbind(setting, t(results))
}


# 37. new fun for rerun in 8.10
get_test_results_new <- function(folder_path, ncores) {
  files <- dir(folder_path)
  map_dfr(files, ~ read_rds(paste0(folder_path, "/", .x)) %>% run_tests_new(ncores))
}
