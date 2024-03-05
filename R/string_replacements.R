var_replacement <- function() {
  c(
    bm_total = "Total biomass",
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
    cpe_int = "Comp effect (Interaction)",
    cpe_env = "Comp effect (Environment)",
    sae_total = "Portfolio effect",
    sae_even = "Portfolio effect (Even)",
    evenness_sae = "Portfolio effect (Evenness)",
    log1_avg_max_int_alive = "Avg interaction strength",
    avg_max_int_alive = "Avg interaction strength",
    avg_int_strength = "Avg interaction strength",
    sd_max_int_alive = "SD interaction strength",
    w_avg_tlvl_alive = "Avg trophic level",
    ct_alive = "Connectance",
    max_tlvl = "Max trophic level",
    w_avg_tlvl = "Avg trophic level",
    avg_omnivory = "Avg omnivory",
    log_richness = "Richness",
    richness = "Richness",
    rho = "Env Corr",
    productivityenrichment = "Enrichment",
    env_stoch = "Env stochasticity",
    Z = "PPMR",
    h = "Hill exponent",
    c = "Predator interference",
    B0 = "Half saturation consumption rate",
    alpha_ij = "Interspecific competition coefficient",
    prod_mass = "Primary producer mass",
    K = "Carrying capacity",
    r = "Maximum growth rate",
    S = "Initial species richness",
    ct = "Initial connectance",
    resp_div = "Response diversity"
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

surround_latex_math <- function(x = "a") {
  paste0("$", x, "$")
}

paste_named_vector <- function(x = NULL, sep = " = ", collapse = ", ") {
  paste0(names(x), sep, x,  collapse = collapse)
}

latex_symbol_replacement <- function() {
  c(
    alpha_ij = "\\alpha_{ij}",
    B0 = "B_0",
    d0 = "d_0",
    r = "r",
    c = "c",
    K = "K",
    eij = "e_{ij}",
    ay = "y",
    ax = "a_x",
    allo_slope = "b",
    prod_mass = "M_p"
  )
}

sim_type_replacement <- function() {
  c(
    "no_rerun" = "Default (5000 timesteps)",
    "rerun" = "Until no extinction",
    "non_allo_no_rerun" = "Fixed death rate",
    "no_rerun_short" = "2000 timesteps",
    "no_rerun_no_disconnected" = "Disconnected species not removed"
  )

}
