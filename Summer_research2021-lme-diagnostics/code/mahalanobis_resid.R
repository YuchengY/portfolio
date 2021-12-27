#' @importFrom diagonals split_vector fatdiag
mahalanobis_ranef.lmerMod <- function(object){
  ngrps <- 0
  mats <- .lmerMod_matrices(object)
  
  n_lev <- length(lme4::getME(object, "flist"))
  
  if(n_lev == 1) {
    Z <- lme4::getME(object, "Z")
    vc <- lme4::VarCorr(object)
    D  <- kronecker(Diagonal(mats$ngrps), bdiag(vc))
    
    eblup <- tcrossprod(D, Z) %*% mats$Vinv %*% resid_marginal(object)
    # vcov_eblup <- D - tcrossprod(D, Z) %*% mats$P %*% tcrossprod(Z, D)
    vcov_eblup <- tcrossprod(D, Z) %*% mats$P %*% tcrossprod(Z, D)
    
    eblup_lst      <- diagonals::split_vector(eblup, size = 2)
    vcov_eblup_lst <- diagonals::fatdiag(vcov_eblup, steps = ngrps(object)) %>%
      diagonals::split_vector(size = 4) %>%
      map(~matrix(.x, nrow = 2, byrow = TRUE))
    
    mah_dist_eblup <- purrr::map2_dbl(eblup_lst, vcov_eblup_lst, ~t(.x) %*% MASS::ginv(.y) %*% .x)
  } else{
    ### Need to check for higher-level models.... D will fail...
  } 
  return(mah_dist_eblup)
}