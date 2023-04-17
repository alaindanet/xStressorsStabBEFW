statistical_averaging_effect <- function(mat) {

  var_sp <- apply(mat, 2, var)
  sd_sp <- apply(mat, 2, sd)

  statistical_averaging <- sum(sd_sp) / sqrt(sum(var_sp))
  statistical_averaging_even <- sqrt(ncol(mat))
  eveness <- statistical_averaging / statistical_averaging_even

  c(
    total = statistical_averaging,
    even = statistical_averaging_even,
    eveness = eveness
  )
}

compensatory_effect <- function(mat) {
  cov_mat <- cov(mat)
  var_sp <- diag(cov_mat)
  var_tot <- sum(cov_mat)
  sqrt(sum(var_sp)) / sqrt(var_tot)
}

asynchrony <- function(mat) {

  cov_mat <- cov(mat)
  var_sp <- diag(cov_mat)
  var_tot <- sum(cov_mat)

  sum(sqrt(var_sp)) / sqrt(var_tot)
}

population_stability <- function(mat) {

  bm_sp <- colMeans(mat)
  bm_tot <- sum(bm_sp)

  sd_sp <- apply(mat, 2, sd)

  bm_tot / sum(sd_sp)
}

community_stability <- function(mat) {

  bm_sp <- colMeans(mat)
  bm_tot <- sum(bm_sp)

  cov_mat <- cov(mat)
  var_tot <- sum(cov_mat)

  bm_tot / sqrt(var_tot)
}

compensatory_effect_env <- function(mat, rho = 0) {
  cov_mat <- cov(mat)
  sd_sp <- sqrt(diag(cov_mat))

  # Build variance-covariance matrix with rho as species cross correlation
  env_mat <- matrix(data = NA, ncol = length(sd_sp), nrow = length(sd_sp))
  for (i in seq_along(sd_sp)) {
    env_mat[i, ] <- sd_sp[i] * sd_sp
  }
  diag_env_mat <- diag(env_mat)
  env_mat <- env_mat * rho
  diag(env_mat) <- diag_env_mat

  sqrt(sum(diag(cov_mat))) / sqrt(sum(env_mat))
}

compensatory_effect_int <- function(mat, rho = 0) {
 compensatory_effect(mat) / compensatory_effect_env(mat, rho = rho)
}
