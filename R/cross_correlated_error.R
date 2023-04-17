# From
# Ripa & Ives (2003): https://www.sciencedirect.com/science/article/pii/S0040580903000893
# Vasseur & Fox (2007): http://doi.wiley.com/10.1111/j.1461-0248.2007.01099.x
# Gouhier et al. (2010): https://www.journals.uchicago.edu/doi/full/10.1086/649579
cross_correlated_response <- function(rho = 0, evt_sd = 1, n = 2, time_len = 100, return_tbl = TRUE) {

  # Independant timeseries
  a <- t(sapply(seq_len(n), function(x) rnorm(time_len, 0, evt_sd)))

  cov_mat <- matrix(rho * evt_sd^2, nrow = nrow(a), ncol = nrow(a))
  diag(cov_mat) <- evt_sd^2

  out <- cov_mat %*% a

  if (return_tbl) {

    dfr <- as.data.frame(t(out))  %>%
      mutate(time = seq_len(ncol(a))) %>%
      pivot_longer(-time, names_to = "species", values_to = "error")
    return(dfr)
  }

  return(out)
}
