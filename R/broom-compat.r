# copied from tidyverse::broom for compatibility
tidy <- function(x, ...) {
  ret <- x[c("estimate", "statistic", "p.value", "parameter")]

  # estimate may have multiple values
  if (length(ret$estimate) > 1) {
    names(ret$estimate) <- paste0("estimate", seq_along(ret$estimate))
    ret <- c(ret$estimate, ret)
    ret$estimate <- NULL

    # special case: in a t-test, estimate = estimate1 - estimate2
    if (x$method == "Welch Two Sample t-test") {
      ret <- c(estimate=ret$estimate1 - ret$estimate2, ret)
    }
  }

  # parameter may have multiple values as well, such as oneway.test
  if (length(x$parameter) > 1) {
    ret$parameter <- NULL
    if (is.null(names(x$parameter))) {
      warning("Multiple unnamed parameters in hypothesis test; dropping them")
    } else {
      message("Multiple parameters; naming those columns ",
              paste(make.names(names(x$parameter)), collapse = ", "))
      ret <- append(ret, x$parameter, after = 1)
    }
  }

  ret <- compact(ret)
  if (!is.null(x$conf.int)) {
    ret <- c(ret, conf.low=x$conf.int[1], conf.high=x$conf.int[2])
  }
  if (!is.null(x$method)) {
    ret <- c(ret, method = as.character(x$method))
  }
  if (!is.null(x$alternative)) {
    ret <- c(ret, alternative = as.character(x$alternative))
  }
  unrowname(as.data.frame(ret))
}

unrowname <- function(x) {
  rownames(x) <- NULL
  x
}

compact <- function(x) Filter(Negate(is.null), x)
