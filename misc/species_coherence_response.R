############################################
#  Test environmental correlated response  #
############################################

library(tidyverse)
library(cowplot)

sdt_e <- 1
rho <- 1

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


# With two species
ti <- tibble(
  rho = c(0, .8),
  noise = map(rho, ~cross_correlated_response(rho = .x, evt_sd = 1, n = 2, time_len = 20))
)


p <- ti %>%
  unnest(noise) %>%
  ggplot(aes(x = time, y = error, color = species)) +
  geom_line() +
  labs(x = "Time", y = "Error") +
  facet_grid(rows = vars(rho)) +
  viridis::scale_color_viridis(discrete = TRUE) +
  cowplot::theme_half_open() +
  theme(legend.position = "none")

ggsave(
  filename = here::here("fig", "coherence_two_sp.pdf"),
  p,
  scale = 1.5,
  width = 88,
  height = 88 * 1.1,
  units = "mm"
)

# With four species
ti <- tibble(
  rho = c(0, .8, -.8),
  noise = map(rho,
    ~cross_correlated_response(rho = .x, evt_sd = 1, n = 4, time_len = 20))
)

p <- ti %>%
  unnest(noise) %>%
  ggplot(aes(x = time, y = error, color = species))+
  geom_line() +
  facet_grid(cols = vars(rho)) +
  cowplot::theme_half_open()

save_plot(
  "~/Téléchargements/coherence_four_sp.png",
  p
)

# Test for trait synchrony
ti <- tibble(
  rho = c(0, .3, -.8),
  noise = map(rho, ~cross_correlated_response(rho = .x, evt_sd = 1, n = 2, time_len = 20))
)

#ti %>%

ti2 <- ti %>%
  mutate(
    cwm = map(noise,
      function (x) {
        x %>%
          pivot_wider(names_from = "species", values_from = "error") %>%
          mutate(
            V1 = V1 + 6,
            V2 = V2 + 3
            ) %>%
          pivot_longer(-time, names_to = "species", values_to = "cwm")
      })
  )

ti2 %>%
  unnest(cwm) %>%
  ggplot(aes(x = time, y = cwm, color = species))+
  geom_line() +
  facet_grid(cols = vars(rho)) +
  cowplot::theme_half_open()

tp <- ti2 %>%
  filter(rho == 0.3) %>%
  unnest(cwm) %>%
  ggplot(aes(x = time, y = cwm, color = species))+
  geom_line() +
  labs(x = "Time", y = "CWM") +
  cowplot::theme_half_open() +
  theme(
    legend.position = "none",
    axis.ticks = element_blank(),
    axis.text = element_blank()
  )
save_plot(
  filename = here::here("fig/trait_synchrony.pdf"),
  tp
)

set.seed(1)
a <- t(sapply(c(3, 6), function(x) rnorm(10, x, 1)))
dfr <- as.data.frame(t(a))  %>%
  mutate(species = seq_len(ncol(a))) %>%
  pivot_longer(-species, names_to = "com", values_to = "trait") %>%
  arrange(com) %>%
  mutate(
    species = seq_len(n()),
    abun = runif(n(), min = .1, max = 1)) %>%
  group_by(com) %>%
  mutate(rel_abun = abun / sum(abun))

library(scales) # ! important
dfr %>%
  ggplot(aes(x = trait, y = rel_abun, fill = com)) +
  geom_col() +
  geom_density(data = dfr,
    aes(x = trait, y = ..density.. / 1.8, color = com, weight = rel_abun),
    adjust = 1, inherit.aes = FALSE)

p_trait <- dfr %>%
  ggplot(aes(x = trait, y = rel_abun, fill = com)) +
  geom_col() +
  geom_density(data = dfr,
    aes(x = trait, y = ..density.. / 2.5, color = com, weight = rel_abun),
    adjust = .8, inherit.aes = FALSE) +
  geom_vline(xintercept = summarise(dfr, avg = sum(rel_abun * trait))$avg) +
  labs(x = "Trait", y = "Relative abundance") +
  cowplot::theme_half_open() +
  #coord_cartesian(expand = FALSE) +
  theme(
    legend.position = "none",
    axis.ticks = element_blank(),
    axis.text = element_blank()
  )

save_plot(
  filename = here::here("fig/p_trait.pdf"),
  p_trait
)
