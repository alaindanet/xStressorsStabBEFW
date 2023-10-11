formula_stab_fw_model <- function(
  resp = "stab_com",
  rand_effect = TRUE,
  term_to_del = "richness : rho : Z +"
  ) {

  fixed_part <- "richness + rho + env_stoch +
  ct_alive + w_avg_tlvl + Z +
  avg_int_strength +
  h +
  richness : rho +
  richness : env_stoch +
  richness : ct_alive +
  richness : avg_int_strength +
  richness : h +
  richness : w_avg_tlvl +
  richness : Z +
  richness : rho : ct_alive +
  richness : rho : w_avg_tlvl +
  richness : rho : avg_int_strength +
  richness : rho : h +
  richness : rho : Z +
  richness : env_stoch : ct_alive +
  richness : env_stoch : Z +
  richness : env_stoch : avg_int_strength +
  richness : env_stoch : h +
  richness : env_stoch : w_avg_tlvl"

  if (rand_effect) {
   rand_part <- "+\n (1|fw_id)"
  } else {
   rand_part <- ""
  }

  form <- paste0(resp, " ~\n", fixed_part, rand_part)
  if (!is.null(term_to_del)) {
    form <- str_remove(form, term_to_del)
  }
  as.formula(form)

}

get_slope_change <- function(
  coeff = fix_fw_stab,
  gen_term = "richness",
  coeff_values = NULL
    ) {

  term <- coeff[str_detect(names(coeff), gen_term)]
  term

  # Get unique term varying with gen_term
  co_varying <- unique(unlist(str_split(names(term), ":")))
  co_varying <- co_varying[!str_detect(co_varying, gen_term)]

  if (!all(co_varying %in% names(coeff_values))) {
    coeff_notin <- co_varying[which(!co_varying %in% names(coeff_values))]
    warning(cat("The coefficients:", coeff_notin, "are not controlled for.\n"))
  }

  # Map over coefficients
  tu <- map_dbl(names(term),
    function(x) {
      # Extract control variables present in each coefficients:
      my_pattern <- paste0("\\b", names(coeff_values), "\\b")
      detected_ctrl_var <- coeff_values[str_detect(x, my_pattern)]

      if (length(detected_ctrl_var) == 0) {
        # if no match, return the coefficient
        out <- term[x]
      } else {
        # if match, multiply the coefficients by the control variables
        out <- reduce(c(term[x], detected_ctrl_var), `*`)
      }
      return(out)
    })

  # Finally sum-up to have the slope values for gen term according to the
  # control variables
  reduce(tu, `+`)
}

model2tibble <- function(x, response = NULL) {
  out <- parameters(x) %>%
    filter(Effects == "fixed", Parameter != "(Intercept)") %>%
    as_tibble() %>%
    janitor::clean_names() %>%
    mutate(
      response = var_replacement()[as.character(formula(x)[2])],
      parameter = str_replace_all(parameter, var_replacement()),
      count_int = str_count(parameter, ":"),
      effect_type = case_when(
        count_int == 0 ~ "single",
        count_int == 1 ~ "double",
        count_int == 2 ~ "triple",
        TRUE ~ "other"

      )
    ) %>%
    select(response, parameter:ci_high, effect_type)
    out
}
plot_effect_type <- function(x) {
  effect <- c("single", "double", "triple")
  p <- map(effect,
    function(y) {
      x %>%
        filter(effect_type == y) %>%
        arrange(coefficient) %>%
        ggplot(aes(
            y = reorder(parameter, -coefficient),
            x = coefficient,xmin = ci_low, xmax = ci_high,
            color = response
            )) +
        geom_vline(xintercept = 0, color = "black", linetype = "dashed", linewidth = 1.5) +
        geom_pointrange(position = position_dodge(width = 0.5)) +
        labs(x = "Standardized Coefficients")
    }
  )
  names(p) <- effect
  p
}

