ggsave_multiple <- function(fns, ...) {
  map(fns, function(x) ggsave(x, ...))
}

rho_cat_colors <- function() {
  c(
    "Low" = "#FF0000",
    "High" = "#02A5E0"
  )
}

resp_div_cat_colors <- function() {
  rho_cat_colors()
}
