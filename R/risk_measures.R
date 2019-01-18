# These functions calculate the two risk measures included in the package, VaR and ES.

# Value at risk ====

create_VaR_func <- function(x, trunc, fitted_pdf, fitted_cdf) {
  function(p = 0.05) {
    if (p < 0 | p > 1)
      stop("p must be between 0 and 1.")
    obj_function <- function(X) fitted_cdf(X) - (1-p)

    # Univariate root finding should be easy: obj_function is a
    # monotonically increasing continuous function. The search interval
    # is initially from the truncation point to the maximum is
    # automatically extended to the right until a room is found.

    lower_limit <- if (is.null(trunc)) 0 else trunc

    uniroot(obj_function, c(lower_limit, max(x)), extendInt = "upX")$root
  }
}

# Expected shortfall ====

create_ES_func <- function(x, trunc, fitted_pdf, fitted_cdf, VaR) {
  function(p = 0.05) {
    if (p < 0 | p > 1)
      stop("p must be between 0 and 1.")
    q0 <- VaR(p)
    fx <- function (x) fitted_pdf(x)*x

    pracma::quadinf(fx, q0, Inf)$Q / pracma::quadinf(fitted_pdf, q0, Inf)$Q
  }
}
