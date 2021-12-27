resid_ranef.lmerMod <- function(object, whichel, standardize = FALSE){
  if(!is.logical(standardize)) {
    stop("standardize must be logical (TRUE or FALSE).")
  }
  
  flist <- lme4::getME(object, "flist")
  
  if(standardize){
    re <- lme4::ranef(object, condVar = TRUE)
    vc <- lme4::VarCorr(object)
    for(i in names(flist)) {
      diag_var <- diag(as.matrix(vc[[i]])) - diag(as.matrix(attr(re[[i]], "postVar")[,,1]))
      re[[i]] <- sweep(re[[i]], 2, sqrt(diag_var), FUN = "/")
    }
  } else {
    re <- lme4::ranef(object)
  }
  
  re[[whichel]]
}


# Simultaneous bands of Schuetzenmesiter and Piepho,
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
                row_names = str_extract(rownames(boot_re$replicates), "\\d+"))
  
  sim_mats <- vector(mode = "list", length = n_re)
  for(i in seq_len(n_re)) {
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
                      })
  
  result <- lapply(bounds, FUN = function(.y) as_tibble(t(.y$Q)) %>% set_names(nm = c("lower", "upper")))
  result <- map2(result, orig_re, .f = ~bind_cols(.x, pred = sort(.y)))
  names(result) <- colnames(res)[seq_len(n_re)]
  
  bind_rows(result, .id = "re")
}


# Plots of the random effects
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
    facet_wrap(~re, ncol  = ncol)
  }
