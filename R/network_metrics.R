gini <- function(x) {
    n <- length(x)
    xx <- sort(x)
    2 * (sum(seq(1, n) * xx)) / (n * sum(xx)) - 1
}
