var_replacement <- function() {
  c(
    stab_com = "Community stability",
    log_stab_com = "Community stability",
    stab_pop = "Population stability",
    pop_stab = "Population stability",
    log_stab_pop = "Population stability",
    log_avg_cv_sp = "Population variability",
    sync = "Synchrony",
    log_sync = "Synchrony",
    log_async = "Asynchrony",
    async = "Asynchrony",
    cpe = "Compensatory effect",
    cpe_int = "Comp. effect (Interaction)",
    cpe_env = "Comp. effect (Environment)",
    sae_total = "Stat. Avg. effect",
    sae_even = "Stat. Avg. effect (Even)",
    evenness_sae = "Stat. Avg. effect (Evenness)",
    log1_avg_max_int_alive = "Interaction strength",
    avg_max_int_alive = "Avg interaction strength",
    avg_int_strength = "Avg interaction",
    sd_max_int_alive = "SD interaction strength",
    w_avg_tlvl_alive = "Avg trophic level",
    ct_alive = "Connectance",
    max_tlvl = "Max. trophic level",
    w_avg_tlvl = "Avg. trophic level",
    avg_omnivory = "Avg. omnivory",
    log_richness = "Richness",
    richness = "Richness",
    rho = "Env Corr",
    productivityenrichment = "Enrichment",
    env_stoch = "Env variance",
    Z = "PPMR"

  )
}
rho_replacement <- function() {
  c(
    `0` = "Maximum",
    `0.5` = "Medium",
    `1` = "No"
  )
}

get_rev_vec_name_val <- function(x = NULL) {
  y <- names(x)
  names(y) <- x
  return(y)
}
