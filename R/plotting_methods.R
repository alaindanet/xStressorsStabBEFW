ggsave_multiple <- function(fns, ...) {
  map(fns, function(x) ggsave(x, ...))
}
