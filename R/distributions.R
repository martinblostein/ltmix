# These functions are concerned with generating the PDF and CDF for each (possibly truncated)
# distribution in a consistent format, as expected by the other functions in the package.

# Probility density functions ====

lognormal_pdf <- function(x, pars, as.log = FALSE) {
    dlnorm(x, meanlog = pars[1], sdlog = pars[2], log = as.log)
}

gamma_pdf <- function(x, pars, as.log=FALSE) {
    if (!all(pars > 0)) {
        Inf
    } else {
        dgamma(x, shape = pars[1], scale = pars[2], log=as.log)
    }
}

weibull_pdf <- function(x, pars, as.log = FALSE) {
    # dweibull sometimes returns NaN instead of zero
    y <- suppressWarnings(dweibull(x, shape = pars[1], scale = pars[2], log=as.log))
    na_replacement <- if(as.log) -Inf else 0
    ifelse(is.nan(y), na_replacement, y)
}

# Cumulative density functions ====

lognormal_cdf <- function(x, pars) {
    plnorm(x, meanlog = pars[1], sdlog = pars[2])
}

gamma_cdf <- function(x, pars) {
    pgamma(x, shape = pars[1], scale = pars[2])
}

weibull_cdf <- function(x, pars) {
    pweibull(x, shape = pars[1], scale = pars[2])
}

# Selection functions ====
# These functions take density names as strings and return the appropriate PDF/CDF

select_pdf <- function(density_name) {
    switch(density_name,
           "lognormal" = lognormal_pdf,
           "gamma" = gamma_pdf,
           "weibull" = weibull_pdf
    )
}


select_cdf <- function(density_name) {
    switch(density_name,
           "lognormal" = lognormal_cdf,
           "gamma" = gamma_cdf,
           "weibull" = weibull_cdf
    )
}

# Truncation ====
# These functions take a density name as a string and a truncation value and return
# the appropriate left-trunctated PDF/CDF

create_trunc_pdf <- function(density_name, trunc) {
    # creates vector-capable truncated pdf, returns as an R function
    # if trunc is null or <= 0, returns the usual pdf

    if (is.null(trunc) || trunc <= 0)
        select_pdf(density_name) # usual pdf
    else {
        pdf <- select_pdf(density_name)
        cdf <- select_cdf(density_name)

        function(x, pars, as.log=FALSE) { # truncated pdf
            if (as.log) {
                ifelse(x >= trunc,
                       pdf(x, pars, as.log = TRUE) - log(1-cdf(trunc, pars)), -Inf)
            } else {
                ifelse(x >= trunc,
                       pdf(x, pars, as.log = FALSE) / (1-cdf(trunc, pars)), 0)
            }
        }
    }
}

create_trunc_cdf <- function(density_name, trunc) {
    # creates vector-capable truncated cdf, returns as an R function
    # if trunc is null or <= 0, returns the usual cdf

    if (is.null(trunc)) #|| trunc <= 0)
        select_cdf(density_name) # usual cdf
    else {
        cdf <- select_cdf(density_name)

        function(x, pars) { # truncated cdf
            ifelse(x >= trunc,
                   (cdf(x, pars) - cdf(trunc, pars)) / (1-cdf(trunc, pars)), 0)
        }
    }
}
