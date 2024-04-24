ggsave_multiple <- function(fns, ...) {
  map(fns, function(x) ggsave(x, ...))
}

rho_cat_colors <- function() {
  c(
    "Low" = "#FF0000",
    "High" = "#02A5E0"
  )
}

stab_color <- function() {
  c(
    "Community stability" = "black", #"#ffffbf"
    "Asynchrony" = "#44AA99",
    "Population stability" = "#E1BE6A"
  )
}

stab_alpha <- function () {
  c(
    "Community stability" = 1,
    "Asynchrony" = 1,
    "Population stability" = 1
  )
}

sem_tot_predictor <- function () {
  c("resp_div", "avg_int_strength", "w_avg_tlvl", "avg_omnivory", "ct_alive",
    "richness", "env_stoch")
}

resp_div_cat_colors <- function() {
  rho_cat_colors()
}

timespecies_mat_to_long_tibble <- function(x) {
  x %>%
    as_tibble() %>%
    mutate(timestep = seq_len(nrow(.))) %>%
    pivot_longer(
      -timestep,
      names_to = "species",
      values_to = "biomass"
    )
}

plot_timeseries_lg_tbl <- function(x, text_font_size = 8) {
  x %>%
    ggplot(aes(x = timestep, y = biomass, color = species)) +
    geom_line() +
    labs(x = "Timesteps", y = "Biomass") +
    theme_half_open(font_size = text_font_size) +
    theme(legend.position = "none")

}

get_stab_fw_label <- function(
  x = sim_for_plot,
  id = 1,
  var = c(
    "richness", "Z", "ct_alive", "w_avg_tlvl", "avg_int_strength",
    "stab_com", "pop_stab", "async", "sae_total", "cpe"
    ),
  sep = ": "
  ) {
  out <- x %>%
    filter(sim_id == id) %>%
    select(where(is.numeric)) %>%
    pivot_longer(everything()) %>%
    deframe()
  out <- out[var]
  paste0(
    var_replacement()[names(out)], sep,
    map_chr(signif(out, digits = 2), ~format(.x, scientific = FALSE)),
    collapse = "\n"
)
}

get_mat_list <- function(x, var = species) {
  out <- x %>%
    pull({{var}})
  names(out) <- x %>%
    pull(sim_id)
  out
}

# Code to plot network
source("https://raw.githubusercontent.com/alaindanet/fishcom/master/R/codeWeb.R")

plot_sim <- function(x = sim_for_plot,
  id = 6827,
  col_palette = NULL,
  node_size_scale = .5,
  font_size = 10,
  tmax = 200,
  replace0by = 1e-6,
  title = NULL,
  var_label = NULL,
  plot_label = TRUE,
  plot_fw = TRUE
  ) {
  # get the good sim
  sim <- x %>%
    filter(sim_id == id)
  # Get the matrix
  mat <- get_mat_list(sim)[[1]]

  if (!is.null(tmax)) {
    mat <- mat[(nrow(mat) - tmax):nrow(mat), ]
  }

  if (!is.null(replace0by)) {
    mat[mat == 0] <- replace0by
  }


  # Get species color
  if (is.null(col_palette)) {
    col_palette <- paletteer::paletteer_dynamic("cartography::multi.pal", ncol(mat))
    col_palette <- setNames(col_palette, paste0("V", seq(1, ncol(mat))))
  }

  # Make the timeseseries plot
  ts_plot <- timespecies_mat_to_long_tibble(mat) %>%
    plot_timeseries_lg_tbl(., text_font_size = font_size) +
    scale_color_manual(values = col_palette) +
    scale_y_continuous(trans = "log10")

  if (!is.null(title)) {
    ts_plot <- ts_plot +
      ggtitle(title)
  }

  # Add the stability and food-web metrics
  if (plot_label) {
    ## Get labels
    if (!is.null(var_label)) {
      stab_fw_label <- get_stab_fw_label(x = x, id = id, var = var_label)
    } else {
      stab_fw_label <- get_stab_fw_label(x = x, id = id)
    }
    ## Add them to the plot
    ts_plot <- ts_plot +
      coord_cartesian(
        expand = FALSE,
        # Do not cut the annotations below
        clip = "off") +
      # This widens the right margin
      theme(plot.margin = unit(c(1, 12, 1, 1), "lines")) +
      annotate("text", x = nrow(mat) + 2, y = 0,
        label = stab_fw_label,
        #https://stackoverflow.com/a/65077171
        size = font_size / .pt - 1,
        hjust = 0,
        lineheight = .8,
        vjust = 0
      )
  }

  # Get the Food-web
  ## Adjency matrix
  int_mat <- get_mat_list(sim, var = int_strength)[[1]]
  adj_mat <- int_mat > 0
  TL <- NetIndices::TrophInd(t(adj_mat))$TL
  ## Plot the food-web
  p <- function(abun = colMeans(mat)) {
    par(
      mar = c(0, 0, 0, 0),
      mgp = c(0, 0, 0),
      pty = "s"
    )
    PlotWeb(
      TL = TL,
      webTL = t(adj_mat),
      colnode = col_palette,
      abund = abun,
      collink = "grey70",
      scale_abun = node_size_scale
    )
  }

  # Combine timeseries and food-web
  if (plot_fw) {
  ggdraw(ts_plot) +
    draw_plot(p,
      x = 1, y = 1,
      scale = .4,
      hjust = 1, vjust = 1,
      halign = 1,
      valign = 1
    )
  } else {
    ts_plot
  }

}
