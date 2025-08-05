summary_distribution <- function(x = NULL, na.rm = FALSE) {

    quant <- quantile(x, probs = c(0, .05, .25, .5, .75, .95, 1), na.rm = na.rm)
    names(quant) <- c("min", "5%", "1st_quart", "median", "2nd_quart", "95%", "max")

    other_desc <- c(
      mean = mean(x, na.rm = na.rm),
      sd = sd(x, na.rm = na.rm),
      n = length(x),
      n_na = length(x[is.na(x)]),
      frac_na = length(x[is.na(x)]) / length(x)
    )

    output <- c(quant, other_desc)
    return(output)
}

get_summary_df <- function(x = NULL, nsignif = 2) {
  x %>%
    pivot_longer(everything(),
      names_to = "var", values_to = "values") %>%
  group_by(var) %>%
  summarise(summ = list(enframe(summary_distribution(values)))) %>%
  unnest(summ) %>%
  pivot_wider(names_from = "name", values_from = "value") %>%
  group_by(var) %>%
  mutate(
    median = format(signif(median, 2), nsmall = 1),
    `Median (Q1, Q3)` = paste0(median, " (",
      format(signif(`1st_quart`, 2), nsmall = 1), ",",
      format(signif(`2nd_quart`, 2), nsmall = 1),
      ")"),
    `Median (5%, 95%)` = paste0(median, " (",
      format(signif(`5%`, 2), nsmall = 1), ",",
      format(signif(`95%`, 2), nsmall = 1),
      ")"),
    `(Min, Max)` = paste0("(",
      format(signif(min, 2), nsmall = 1, scientific = FALSE),
      ",",
      format(signif(max, 2), nsmall = 1),
      ")")
  )
}

get_sim_summary_stat <- function(x = sim_fw) {
  get_summary_df(
    x = x %>% select(where(is.double))
    ) %>%
    mutate(
      type = map_chr(var,
        function (x) {
          if (x %in% c("stab_com", "pop_stab", "async", "sae_doak", "cpe_int", "cpe_env")) {
            "stability"
          } else if (x %in% c("avg_int_strength", "richness", "ct_alive", "w_avg_tlvl", "max_tlvl", "avg_omnivory")) {
            "food-web"
          } else {
            NA
          }
        })
      ) %>%
  filter(!is.na(type)) %>%
  arrange(desc(type))
}
