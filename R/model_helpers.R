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

