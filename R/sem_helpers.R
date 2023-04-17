get_table_semEff <- function(semeff = NULL) {

  effect_type <-  c("^DIRECT", "^INDIRECT", "TOTAL", "MEDIATORS")
  x  <- semeff[-1, ]
  effect_type_col <- x[[1]]
  predictor <- x[[2]]

  row_effect_type <- map_int(effect_type, ~which(str_detect(x[[1]], .x)))

  numeric_table <- map(x, as.numeric)
  # suppr table columns
  mask_chr_column <- map_lgl(numeric_table, ~all(is.na(.x)))
  # Keep the second columns:
  mask_chr_column[2] <- FALSE

  tab <- x[, !mask_chr_column]
  colnames(tab)[1] <- "predictor"

  output <- rbind(
    cbind(effect_type = "direct", tab[row_effect_type[1]:row_effect_type[2]-1, ]),
    cbind(effect_type = "indirect", tab[row_effect_type[2]:row_effect_type[3]-1, ]),
    cbind(effect_type = "total", tab[row_effect_type[3]:row_effect_type[4]-1, ]),
    cbind(effect_type = "mediators", tab[row_effect_type[4]:nrow(tab), ])
  )

  numeric_rows <- str_detect(output[["Effect"]], "\\d")
  output <- output[numeric_rows,]
  output[, -c(1:2)] <- apply(output[, -c(1:2)], 2, as.numeric)

  janitor::clean_names(output)
}

from_semEff_to_table <- function(x = NULL) {

  output <- map_dfr(x$Summary[-1], get_table_semEff, .id = "response") %>%
    as_tibble()

  output[, c("response", "predictor")] <- apply(
    output[, c("response", "predictor")], 2,
    function(x) str_replace_all(x, c("\\." = "_", "\\s" = ""))
  )
  output[, colnames(output) != "bias"]

}
