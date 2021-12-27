# Matrix sqrt needed for lesaffre-verbeke stat in Singer et al.
matrix_sqrt <- function(mat) {
  mat <- as.matrix(mat)
  decomp <- eigen(mat)
  V <- decomp$vectors
  S <- sqrt(decomp$values)
  V %*% diag(S) %*% t(V)
}

# Currently only for 2-level models fit via lme4
lesaffre_verbeke <- function(object) {
  mats <- HLMdiag:::.lmerMod_matrices(object)
  mar_var <- mats$V - mats$X %*% tcrossprod(mats$XVXinv, mats$X)
  mar_resid <- HLMdiag:::resid_marginal(object, type = "raw")
  
  grp_index <- getME(object, "flist")[[1]]
  grp_names <- levels(grp_index)
  grp_sizes <- table(grp_index)
  
  diagnostic <- vector("numeric", length(grp_names))
  for(i in seq_along(grp_names)) {
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
