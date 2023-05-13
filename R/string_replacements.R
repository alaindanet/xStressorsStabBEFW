var_replacement <- function() {
  c(
    stab_com = "Community stability",
    log_stab_com = "Community stability",
    stab_pop = "Population stability",
    log_stab_pop = "Population stability",
    log_avg_cv_sp = "Population variability",
    sync = "Synchrony",
    log_sync = "Synchrony",
    log_async = "Asynchrony",
    async = "Asynchrony",
    log1_avg_max_int_alive = "Interaction strength",
    ct_alive = "Connectance",
    log_richness = "Richness",
    rho = "Env Corr",
    productivityenrichment = "Enrichment",
    env_stoch = "Env variance"

  )
}
rho_replacement <- function() {
  c(
    `0` = "Maximum",
    `0.5` = "Medium",
    `1` = "No"
  )

}
