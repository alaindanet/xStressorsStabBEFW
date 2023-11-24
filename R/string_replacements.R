var_replacement <- function() {
  c(
    stab_com = "Community stability",
    log_stab_com = "Community stability",
    stab_pop = "Population stability",
    pop_stab = "Population stability",
    log_stab_pop = "Population stability",
    log_avg_cv_sp = "Population variability",
    async = "Asynchrony",
    log_async = "Asynchrony",
    sync = "Synchrony",
    log_sync = "Synchrony",
    cpe = "Compensatory effect",
    cpe_int = "Comp effect (Interaction)",
    cpe_env = "Comp effect (Environment)",
    sae_total = "Stat Avg effect",
    sae_even = "Stat Avg effect (Even)",
    evenness_sae = "Stat Avg effect (Evenness)",
    richness = "Richness",
    log_richness = "Richness",
    avg_int_strength = "Avg interaction strength",
    avg_max_int_alive = "Avg interaction strength",
    log1_avg_max_int_alive = "Avg interaction strength",
    sd_max_int_alive = "SD interaction strength",
    w_avg_tlvl = "Avg trophic level",
    w_avg_tlvl_alive = "Avg trophic level",
    ct_alive = "Connectance",
    max_tlvl = "Max trophic level",
    avg_omnivory = "Avg omnivory",
    productivityenrichment = "Enrichment",
    resp_div = "Response diversity",
    rho = "Env Corr",
    env_stoch = "Env variance",
    Z = "PPMR",
    h = "Hill exponent",
    S = "Initial species richness",
    ct = "Initial connectance"
  )
}

get_term_replacement <- function() {
  c(":" = " x\n")
}

rho_replacement <- function() {
  c(
    `0` = "High",
    `0.5` = "Medium",
    `0.75` = "Low",
    `1` = "No"
  )
}

resp_div_replacement <- function() {
  c(
    `1` = "High",
    `0.5` = "Medium",
    `0.25` = "Low",
    `0` = "No"
  )
}

get_rev_vec_name_val <- function(x = NULL) {
  y <- names(x)
  names(y) <- x
  return(y)
}
