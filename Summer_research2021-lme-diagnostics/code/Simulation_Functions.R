# This R Script is for simulation functions (e.g sim.hlmvar(),sim.hlmlin()) used in the project


# 1. sim.hlmlin
# 2. sim.hlmvar
# 3. sim.hlm
# 4. sim.hlm_re
# 5. sim.hlmfix
# 6. sim.hlm_heterlin
# 7. sim.hlmvar_norm
# 8. sim.hlm_linnorm
# 9. sim.hlm_spec


# 1. Simulate models with linear / nonlinear (quadratic) fixed effect terms
sim.hlmlin <- function(simdata.df, rel.dsn, var) {
  sim_design <- simdata.df # simdata.df Contain id, grouping and two predictors

  b0 <- 5 # Intercept
  b1 <- 1 # X_ij slope
  b2 <- 1 # Z_ij slope

  if (var == "large_error") {
    sd_error <- 5 # large error
    sd_int <- 1 # small random_effect
    sd_slope <- 1
  }

  if (var == "large_re") {
    sd_error <- 1 # small error
    sd_int <- 5 # large random_effect
    sd_slope <- 5
  }

  rho <- 0.5
  cov.matrix <- matrix(c(sd_int^2, sd_int * sd_slope * rho, sd_int * sd_slope * rho, sd_slope^2),
    nrow = 2,
    byrow = TRUE
  )

  m <- max(sim_design$group) # number of groups
  re_df <- data.frame(cbind(group = 1:m, mvrnorm(m, c(0, 0), Sigma = cov.matrix)))
  colnames(re_df) <- c("group", "intercept", "slope")

  # Linearity
  if (rel.dsn == "linear") { # no change
    sim_df <- left_join(sim_design, re_df, by = "group") %>%
      mutate(
        error = rnorm(nrow(.), mean = 0, sd = sd_error),
        Y_ij = b0 + b1 * X_ij + b2 * Z_ij + intercept + slope * X_ij + error
      )
  }

  if (rel.dsn == "sq") { # squaring Z_ij
    sim_df <- left_join(sim_design, re_df, by = "group") %>%
      mutate(
        error = rnorm(nrow(.), mean = 0, sd = sd_error),
        Y_ij = b0 + b1 * X_ij + b2 * (Z_ij)^2 + intercept + slope * X_ij + error
      )
  }

  # y.df <- as.data.frame( as.matrix(sim_df$Y_ij) )
  # colnames(y.df) <- paste("sim_", 1:ncol(y.df), sep = "")
  return(sim_df)
}



# 2. Simulate different levels of non constant variance / heteroscedasticity
sim.hlmvar <- function(simdata.df, evar.dsn, var) {
  sim_design <- simdata.df # simdata.df Contain id, grouping and two predictors

  b0 <- 5 # Intercept
  b1 <- 1 # X_ij slope
  b2 <- 1 # Z_ij slope

  if (var == "large_error") {
    sd_error <- 5 # large error
    sd_int <- 1 # small random_effect
    sd_slope <- 1
  }

  if (var == "large_re") {
    sd_error <- 1 # small error
    sd_int <- 5 # large random_effect
    sd_slope <- 5
  }

  evar.dsn <- as.numeric(evar.dsn)

  rho <- 0.5
  cov.matrix <- matrix(c(sd_int^2, sd_int * sd_slope * rho, sd_int * sd_slope * rho, sd_slope^2),
    nrow = 2,
    byrow = TRUE
  )

  m <- max(sim_design$group) # number of groups
  re_df <- data.frame(cbind(group = 1:m, mvrnorm(m, c(0, 0), Sigma = cov.matrix)))
  colnames(re_df) <- c("group", "intercept", "slope")

  sim_df <- left_join(sim_design, re_df, by = "group") %>%
    mutate(
      error = rnorm(nrow(.), mean = 0, sd = sqrt(sd_error^2 + evar.dsn * (X_ij - min(X_ij)))),
      Y_ij = b0 + b1 * X_ij + b2 * Z_ij + intercept + slope * X_ij + error
    )

  # y.df <- as.data.frame( as.matrix(sim_df$Y_ij) )
  # colnames(y.df) <- paste("sim_", 1:ncol(y.df), sep = "")

  return(sim_df)
}



# 3. sim.hlm
sim.hlm <- function(simdata.df, e.dsn, var) {
  sim_design <- simdata.df # Contain id, grouping and two predictors

  sim_design <- simdata.df # simdata.df Contain id, grouping and two predictors

  b0 <- 5 # Intercept
  b1 <- 1 # X_ij slope
  b2 <- 1 # Z_ij slope

  if (var == "large_error") {
    sd_error <- 5 # large error
    sd_int <- 1 # small random_effect
    sd_slope <- 1
  }

  if (var == "large_re") {
    sd_error <- 1 # small error
    sd_int <- 5 # large random_effect
    sd_slope <- 5
  }

  rho <- 0.5

  m <- max(sim_design$group) # number of groups
  n_obs <- nrow(sim_design) # number of observations

  # cov matrix
  cov.matrix <- matrix(c(sd_int^2, sd_int * sd_slope * rho, sd_int * sd_slope * rho, sd_slope^2),
    nrow = 2,
    byrow = TRUE
  )

  # random effect
  re_df <- data.frame(cbind(group = 1:m, mvrnorm(m, c(0, 0), Sigma = cov.matrix)))
  colnames(re_df) <- c("group", "intercept", "slope")

  # skewed intercept
  if (e.dsn == "norm") {
    e <- rnorm(n = n_obs, mean = 0, sd = sd_error)
    # e <- matrix(e, nc = n_obs)
    sim_df <- left_join(data.frame(sim_design), re_df, by = "group") %>%
      mutate(Y_ij = b0 + b1 * X_ij + b2 * Z_ij + intercept + slope * X_ij + e)
  }


  if (e.dsn == "skewness_3") {
    # fact: kurtosis>1+skewness^2

    # extremely skewed: skewness = 3, kurtosis = 11
    p <- c(mean = 0, variance = sd_error^2, skewness = 3, kurtosis = 11)
    e <- rpearson(n_obs, moments = p)
    # e <- matrix(e, nc = n_obs)
    sim_df <- left_join(data.frame(sim_design), re_df, by = "group") %>%
      mutate(Y_ij = b0 + b1 * X_ij + b2 * Z_ij + intercept + slope * X_ij + e)
  }

  # moderately skewed: skewness = 1.5, kurtosis = 6
  if (e.dsn == "skewness_1.5") {
    p <- c(mean = 0, variance = sd_error^2, skewness = 1.5, kurtosis = 6)
    e <- rpearson(n_obs, moments = p)
    # e <- matrix(e, nc = n_obs)
    sim_df <- left_join(data.frame(sim_design), re_df, by = "group") %>%
      mutate(Y_ij = b0 + b1 * X_ij + b2 * Z_ij + intercept + slope * X_ij + e)
  }

  # slightlyy skewed: skewness = 0.8, kurtosis = 4
  if (e.dsn == "skewness_0.8") {
    p <- c(mean = 0, variance = sd_error^2, skewness = 0.8, kurtosis = 4)
    e <- rpearson(n_obs, moments = p)
    # e <- matrix(e, nc = n_obs)
    sim_df <- left_join(data.frame(sim_design), re_df, by = "group") %>%
      mutate(Y_ij = b0 + b1 * X_ij + b2 * Z_ij + intercept + slope * X_ij + e)
  }


  if (e.dsn == "bimodal") {
    mu1 <- -1 # so their mean cancels out to 0
    mu2 <- 1
    sig1 <- sd_error
    sig2 <- sd_error
    cpct <- 0.4 # take 40% from one and 60% from another
    e <- bimodalDistFunc(n_obs, cpct, mu1, mu2, sig1, sig2)
    # e <- matrix(e, nc = n_obs)
    sim_df <- left_join(data.frame(sim_design), re_df, by = "group") %>%
      mutate(Y_ij = b0 + b1 * X_ij + b2 * Z_ij + intercept + slope * X_ij + e)
  }

  return(sim_df)
}


# 4. sim.hlm_re
sim.hlm_re <- function(simdata.df, re.dsn, var) {
  sim_design <- simdata.df # Contain id, grouping and two predictors

  sim_design <- simdata.df # simdata.df Contain id, grouping and two predictors

  b0 <- 5 # Intercept
  b1 <- 1 # X_ij slope
  b2 <- 1 # Z_ij slope

  if (var == "large_error") {
    sd_error <- 5 # large error
    sd_int <- 1 # small random_effect
    sd_slope <- 1
  }

  if (var == "large_re") {
    sd_error <- 1 # small error
    sd_int <- 5 # large random_effect
    sd_slope <- 5
  }

  rho <- 0.5

  # n<-nrow(data.frame(sim_design)) # number of obs
  m <- max(sim_design$group) # number of groups

  # cov matrix for normal settings
  cov.matrix <- matrix(c(sd_int^2, sd_int * sd_slope * rho, sd_int * sd_slope * rho, sd_slope^2),
    nrow = 2,
    byrow = TRUE
  )

  # norm
  if (re.dsn == "norm_re") {
    re_df <- data.frame(cbind(group = 1:m, mvrnorm(m, c(0, 0), Sigma = cov.matrix)))
    colnames(re_df) <- c("group", "intercept", "slope")

    sim_df <- left_join(data.frame(sim_design), re_df, by = "group") %>%
      mutate(
        error = rnorm(nrow(.), mean = 0, sd = sd_error),
        Y_ij = b0 + b1 * X_ij + b2 * Z_ij + intercept + slope * X_ij + error
      )
    y.df <- as.data.frame(as.matrix(sim_df$Y_ij))
    colnames(y.df) <- paste("sim_", 1:ncol(y.df), sep = "")
  }


  # skewed intercept
  if (re.dsn == "mildly_skewed_re_intercept") {
    re_df <- data.frame(cbind(group = 1:m, bvn_intercept(m, 0, sd_int, 0, sd_slope, rho)))
    colnames(re_df) <- c("group", "intercept", "slope")

    sim_df <- left_join(data.frame(sim_design), re_df, by = "group") %>%
      mutate(
        error = rnorm(nrow(.), mean = 0, sd = sd_error),
        Y_ij = b0 + b1 * X_ij + b2 * Z_ij + intercept + slope * X_ij + error
      )
    y.df <- as.data.frame(as.matrix(sim_df$Y_ij))
    colnames(y.df) <- paste("sim_", 1:ncol(y.df), sep = "")
    # return( y.df )
  }

  # skewed slope
  if (re.dsn == "mildly_skewed_re_slope") {
    re_df <- data.frame(cbind(group = 1:m, bvn_slope(m, 0, sd_int, 0, sd_slope, rho)))
    colnames(re_df) <- c("group", "intercept", "slope")

    sim_df <- left_join(data.frame(sim_design), re_df, by = "group") %>%
      mutate(
        error = rnorm(nrow(.), mean = 0, sd = sd_error),
        Y_ij = b0 + b1 * X_ij + b2 * Z_ij + intercept + slope * X_ij + error
      )
  }
  return(sim_df)
}


# 5. simulating models with one fixed effect omitted
sim.hlmfix <- function(simdata.df, z.dsn, var) {
  sim_design <- simdata.df # Contain id, grouping and two predictors

  sim_design <- simdata.df # simdata.df Contain id, grouping and two predictors

  b0 <- 5 # Intercept
  b1 <- 1 # X_ij slope
  b2 <- 1 # Z_ij slope

  if (var == "large_error") {
    sd_error <- 5 # large error
    sd_int <- 1 # small random_effect
    sd_slope <- 1
  }

  if (var == "large_re") {
    sd_error <- 1 # small error
    sd_int <- 5 # large random_effect
    sd_slope <- 5
  }

  rho <- 0.5
  cov.matrix <- matrix(c(sd_int^2, sd_int * sd_slope * rho, sd_int * sd_slope * rho, sd_slope^2),
    nrow = 2,
    byrow = TRUE
  )

  m <- max(sim_design$group) # number of groups
  re_df <- data.frame(cbind(group = 1:m, mvrnorm(m, c(0, 0), Sigma = cov.matrix)))
  colnames(re_df) <- c("group", "intercept", "slope")

  if (z.dsn == "full") { 
    sim_df <- left_join(sim_design, re_df, by = "group") %>%
      mutate(
        error = rnorm(nrow(.), mean = 0, sd = sd_error),
        Y_ij = b0 + b1 * X_ij + b2 * Z_ij + intercept + slope * X_ij + error
      )
  }

  if (z.dsn == "reduced") { 
    sim_df <- left_join(sim_design, re_df, by = "group") %>%
      mutate(
        error = rnorm(nrow(.), mean = 0, sd = sd_error),
        Y_ij = b0 + b1 * X_ij +  b2 * Z_ij + intercept + slope * X_ij + error
      )
  }

  # y.df <- as.data.frame( as.matrix(sim_df$Y_ij) )
  # colnames(y.df) <- paste("sim_", 1:ncol(y.df), sep = "")

  return(sim_df)
}



# 6. simulating models with one fixed effect squared and non constant variance
sim.hlm_heterlin <- function(simdata.df, helin.dsn, var) {
  sim_design <- simdata.df # simdata.df Contain id, grouping and two predictors

  b0 <- 5 # Intercept
  b1 <- 1 # X_ij slope
  b2 <- 1 # Z_ij slope

  if (var == "large_error") {
    sd_error <- 5 # large error
    sd_int <- 1 # small random_effect
    sd_slope <- 1
  }

  if (var == "large_re") {
    sd_error <- 1 # small error
    sd_int <- 5 # large random_effect
    sd_slope <- 5
  }

  rho <- 0.5
  cov.matrix <- matrix(c(sd_int^2, sd_int * sd_slope * rho, sd_int * sd_slope * rho, sd_slope^2),
    nrow = 2,
    byrow = TRUE
  )

  m <- max(sim_design$group) # number of groups
  re_df <- data.frame(cbind(group = 1:m, mvrnorm(m, c(0, 0), Sigma = cov.matrix)))
  colnames(re_df) <- c("group", "intercept", "slope")

  if (helin.dsn == "linear_homo") { # no change
    sim_df <- left_join(sim_design, re_df, by = "group") %>%
      mutate(
        error = rnorm(nrow(.), mean = 0, sd = sqrt(sd_error^2 + 0 * (X_ij - min(X_ij)))),
        Y_ij = b0 + b1 * X_ij + b2 * Z_ij + intercept + slope * X_ij + error
      )
  }

  if (helin.dsn == "sq_2") {
    sim_df <- left_join(sim_design, re_df, by = "group") %>%
      mutate(
        error = rnorm(nrow(.), mean = 0, sd = sqrt(sd_error^2 + 2 * (X_ij - min(X_ij)))),
        Y_ij = b0 + b1 * X_ij + b2 * (Z_ij)^2 + intercept + slope * X_ij + error
      )
  }

  if (helin.dsn == "sq_4") {
    sim_df <- left_join(sim_design, re_df, by = "group") %>%
      mutate(
        error = rnorm(nrow(.), mean = 0, sd = sqrt(sd_error^2 + 4 * (X_ij - min(X_ij)))),
        Y_ij = b0 + b1 * X_ij + b2 * (Z_ij)^2 + intercept + slope * X_ij + error
      )
  }

  if (helin.dsn == "sq_8") {
    sim_df <- left_join(sim_design, re_df, by = "group") %>%
      mutate(
        error = rnorm(nrow(.), mean = 0, sd = sqrt(sd_error^2 + 8 * (X_ij - min(X_ij)))),
        Y_ij = b0 + b1 * X_ij + b2 * (Z_ij)^2 + intercept + slope * X_ij + error
      )
  }

  # y.df <- as.data.frame( as.matrix(sim_df$Y_ij) )
  # colnames(y.df) <- paste("sim_", 1:ncol(y.df), sep = "")
  return(sim_df)
}


# 7. simulating models with skewed error term and non-constant variance
sim.hlmvar_norm <- function(simdata.df, evar.dsn, var) {
  sim_design <- simdata.df # Contain id, grouping and two predictors
  b0 <- 5 # Intercept
  b1 <- 1 # X_ij slope
  b2 <- 1 # Z_ij slope

  if (var == "large_error") {
    sd_error <- 5 # large error
    sd_int <- 1 # small random_effect
    sd_slope <- 1
  }

  if (var == "large_re") {
    sd_error <- 1 # small error
    sd_int <- 5 # large random_effect
    sd_slope <- 5
  }

  rho <- 0.5

  m <- max(sim_design$group) # number of groups

  # cov matrix for random effects
  cov.matrix <- matrix(c(sd_int^2, sd_int * sd_slope * rho, sd_int * sd_slope * rho, sd_slope^2),
    nrow = 2,
    byrow = TRUE
  )

  # generate random effects
  re_df <- data.frame(cbind(group = 1:m, mvrnorm(m, c(0, 0), Sigma = cov.matrix)))
  colnames(re_df) <- c("group", "intercept", "slope")


  # skewed + hetero

  # mildly skewed
  p_df <- c(mean = 0, variance = sd_error^2 + parse_number(evar.dsn) * (sim_design$X_ij - min(sim_design$X_ij)), skewness = 1.5, kurtosis = 6)
  p <- vec_p2(m, p_df)

  error <- vec.rpearson(1, moments = p)

  sim_df <- left_join(sim_design, re_df, by = "group") %>%
    mutate(Y_ij = b0 + b1 * X_ij + b2 * Z_ij + intercept + slope * X_ij + error)

  return(sim_df)
}



# 9. sim function for non-linear and non-normal
sim.hlm_linnorm <- function(simdata.df, lin_norm.dsn, var) {
  sim_design <- simdata.df # simdata.df Contain id, grouping and two predictors

  b0 <- 5 # Intercept
  b1 <- 1 # X_ij slope
  b2 <- 1 # Z_ij slope

  if (var == "large_error") {
    sd_error <- 5 # large error
    sd_int <- 1 # small random_effect
    sd_slope <- 1
  }

  if (var == "large_re") {
    sd_error <- 1 # small error
    sd_int <- 5 # large random_effect
    sd_slope <- 5
  }

  rho <- 0.5 # moderate correlation
  cov.matrix <- matrix(c(sd_int^2, sd_int * sd_slope * rho, sd_int * sd_slope * rho, sd_slope^2),
    nrow = 2,
    byrow = TRUE
  )

  m <- max(sim_design$group) # number of groups
  n_obs <- nrow(sim_design) # number of observations


  re_df <- data.frame(cbind(group = 1:m, mvrnorm(m, c(0, 0), Sigma = cov.matrix)))
  colnames(re_df) <- c("group", "intercept", "slope")


  if (lin_norm.dsn == "reduced_skew") {
    p <- c(mean = 0, variance = sd_error^2, skewness = 1.5, kurtosis = 6)
    e <- rpearson(n_obs, moments = p)


    sim_df <- left_join(sim_design, re_df, by = "group") %>%
      mutate(
        Y_ij = b0 + b1 * X_ij + b2 * (Z_ij)^2 + intercept + slope * X_ij + e
      )
  }

  if (lin_norm.dsn == "reduced_bimodal") {
    mu1 <- -1 # so their mean cancels out to 0
    mu2 <- 1
    sig1 <- sd_error
    sig2 <- sd_error
    cpct <- 0.4 # take 40% from one and 60% from another

    e <- bimodalDistFunc(n_obs, cpct, mu1, mu2, sig1, sig2)

    sim_df <- left_join(sim_design, re_df, by = "group") %>%
      mutate(Y_ij = b0 + b1 * X_ij + b2 * (Z_ij)^2 + intercept + slope * X_ij + e)
  }

  return(sim_df)
}


# 9. sim function for special cases
sim.hlm_spec <- function(simdata.df, var, special, J) {
  sim_design <- simdata.df # Contain id, grouping and two predictors

  b0 <- 5 # Intercept
  b1 <- 1 # X_ij slope
  b2 <- 1 # Z_ij slope

  if (var == "large_error") {
    sd_error <- 5 # large error
    sd_int <- 1 # small random_effect
    sd_slope <- 1
  }

  if (var == "large_re") {
    sd_error <- 1 # small error
    sd_int <- 5 # large random_effect
    sd_slope <- 5
  }

  rho <- 0.5

  m <- max(sim_design$group) # number of groups
  n_obs <- nrow(sim_design) # number of observations


  if (special == "time_seq" | special == "ar_error") {
    # cov matrix
    cov.matrix <- matrix(c(sd_int^2, sd_int * sd_slope * rho, sd_int * sd_slope * rho, sd_slope^2),
      nrow = 2,
      byrow = TRUE
    )

    # random effect
    re_df <- data.frame(cbind(group = 1:m, mvrnorm(m, c(0, 0), Sigma = cov.matrix)))
    colnames(re_df) <- c("group", "intercept", "slope")
  }


  if (special == "ar_error") {
    ar.val <- 0.4
    # error sd
    # sigma <- 1.5
    # ar error by group
    ### note: use arima.sim(model=list(ar=ar.val), n=x) * sqrt(1-ar.val^2) * sigma
    ### construction, so that the true error SD is equal to sigma

    eij <- unlist(sapply(J, function(x) arima.sim(model = list(order = c(1, 0, 0), ar = ar.val), n = x) * sqrt(1 - ar.val^2) * sd_error))
    # weird thing, if group is balanced, we have a vector; if unbalanced, we get a matrix/array
    filter <- class(eij) == "matrix"
    if (filter[1] == TRUE) {
      e <- concat(m, eij)
    }
    if (filter[1] != TRUE) {
      e <- eij
    }

    sim_df <- left_join(data.frame(sim_design), re_df, by = "group") %>%
      mutate(Y_ij = b0 + b1 * X_ij + b2 * Z_ij + intercept + slope * X_ij + e)
  }

  if (special == "time_seq") {
    b3 <- 0.5
    e <- rnorm(n = n_obs, mean = 0, sd = sd_error)

    sim_df <- left_join(data.frame(sim_design), re_df, by = "group") %>%
      mutate(Y_ij = b0 + b1 * X_ij + b3 * T_ij + b2 * Z_ij + intercept + slope * X_ij + e)
  }

  if (special == "re_int") {
    re_df <- as.data.frame(cbind(group = 1:m, rnorm(m, mean = 0, sd = sd_int)))
    colnames(re_df) <- c("group", "intercept")

    e <- rnorm(n = n_obs, mean = 0, sd = sd_error)
    sim_df <- left_join(as.data.frame(sim_design), re_df, by = "group") %>%
      mutate(Y_ij = b0 + b1 * X_ij + b2 * Z_ij + intercept + e)
  }

  return(sim_df)
}
