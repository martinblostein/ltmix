
#' Fit a Left-truncated mixture model (LTMM)
#'
#' @param x data vector
#' @param G number of components
#' @param distributions densities to combine
#' @param trunc left truncation point (optional)
#' @param EM_init_method  initialization method for EM algorithm
#' @param EM_starts number of random starts for initialization of EM algorithm. (only for G > 1)
#' @param init_pars initial parameter values (list of length G)
#' @param init_pi manually specified initial component proportions (for init_method=specified)
#' @param init_classes manually specified initial classes. will overwrite init_pars and init_pi
#' @param one_group_reps number of random starts for each numerical optimization in 1-component model
#' @param eps stopping tolerance for EM algoithm
#' @param max.it maximum number of iterations of EM algorithm
#' @param verbose print information as fitting progresses?
#' @param parallel fit models in parallel?
#' @param cores number of processes used for parallel computation. if NULL detect.cores() used
#' @param save_each_fitsave each model as it is produced, in a time-stamped directory (safer)
#'
#' @description This function generate a mixture of lognormal, gamma, and weibull distributions
#'
#' @details Blostein Martina and Miljkovic Tatjana (2019). On On modeling left-truncated loss data using mixtures of distributions. Insurance: Mathematics and Economics. DOI:10.1016/j.insmatheco.2018.12.001
#'
#' @return an ltmmCombo model object
#'
#' @export
ltmmCombo <- function(x, G, distributions = c("lognormal", "gamma", "weibull"), trunc = NULL,
                      EM_init_method = "emEM", EM_starts = 5, init_pars = NULL, init_pi = NULL,
                      init_classes = NULL, one_group_reps = 50, eps = 1e-6, max.it = 1000,
                      verbose = FALSE, parallel = FALSE, cores = NULL, save_each_fit = FALSE) {

    combos <- gtools::combinations(length(distributions), G, distributions, repeats.allowed = T)

    if (save_each_fit) {
        results.dir <- paste0("fits-", gsub(" ", "-", gsub("-|:", "", Sys.time())))
        dir.create(results.dir, showWarnings = FALSE)
    }

    if (parallel) {
        verbose <- FALSE

        if (is.null(cores))
            cores <- parallel::detectCores()

        all.fits <- parallel::mclapply(1:nrow(combos), function(i) {
            tryCatch(
                expr = {
                    w <- ltmm(x, G, combos[i,], trunc, EM_init_method, EM_starts, init_pars[[i]],
                              init_pi[[i]], init_classes, one_group_reps, eps, max.it, verbose)
                    if (save_each_fit)
                        saveRDS(w, file.path(results.dir, paste0("fit_", i, ".rds")))
                    w
                },
                error = function(w) {
                    list(bic = Inf, aic = Inf, ll = -Inf, error=w, densities=combos[i,])
                }
            )
        }, mc.cores = cores)
    } else {
        all.fits <- lapply(1:nrow(combos), function(i) {
            tryCatch(
                expr = {
                    w <- ltmm(x, G, combos[i,], trunc, EM_init_method, EM_starts, init_pars[[i]],
                              init_pi[[i]], init_classes, one_group_reps, eps, max.it, verbose)
                    if (save_each_fit)
                        saveRDS(w, file.path(results.dir, paste0("fit_", i, ".rds")))
                    w
                }, error = function(w) {
                    list(bic = Inf, aic = Inf, ll = -Inf, error=w, densities=combos[i,])
                }
            )
        })
    }

    all.bic <- sapply(all.fits, getElement, 'bic')
    best.bic.ind <- which.min(all.bic)

    all.aic <- sapply(all.fits, getElement, 'aic')
    best.aic.ind <- which.min(all.aic)

    all.ll <- sapply(all.fits, getElement, 'll')
    all.VaR <- function(p = 0.05)
        sapply(all.fits, function(fit) fit$VaR(p))
    all.ES <- function(p = 0.05)
        sapply(all.fits, function(fit) fit$ES(p))

    summary_table <- data.frame(bic = all.bic,
                                aic = all.aic,
                                ll = all.ll,
                                Var_95 = tryCatch(all.VaR(0.05), error = function(e) NA),
                                Var_99 = tryCatch(all.VaR(0.01), error = function(e) NA),
                                ES_95  = tryCatch(all.ES(0.05),  error = function(e) NA),
                                ES_99  = tryCatch(all.ES(0.01),  error = function(e) NA))

    summary_table <- cbind(unname(combos), summary_table)
    names(summary_table)[1:G] = paste0("G", 1:G)

    output <- list(
        x = x,
        G = G,
        distributions = distributions,
        combos = combos,
        all.fits = all.fits,

        all.bic = all.bic,
        best.bic.fit = all.fits[[best.bic.ind]],
        best.bic = all.bic[best.bic.ind],
        best.bic.combo = combos[best.bic.ind,],

        all.aic = all.aic,
        best.aic.fit = all.fits[[best.aic.ind]],
        best.aic = all.aic[best.aic.ind],
        best.aic.combo = combos[best.aic.ind,],

        all.ll = all.ll,

        summary_table = summary_table

    )

    class(output) <- "ltmmCombo"

    if (save_each_fit)
        saveRDS(output, file.path(results.dir, "all_fits.rds"))

    output
}

# Main function for fitting an individual model:
# ltmm = "Left Truncated Mixture Models"
#' Fit a Left-truncated mixture model (LTMM)
#'
#' @param x data vector
#' @param G number of components
#' @param distributions densities to combine
#' @param trunc left truncation point (optional)
#' @param EM_init_method  initialization method for EM algorithm
#' @param EM_starts number of random starts for initialization of EM algorithm. (only for G > 1)
#' @param init_pars initial parameter values (list of length G)
#' @param init_pi manually specified initial component proportions (for init_method=specified)
#' @param init_classes manually specified initial classes. will overwrite init_pars and init_pi
#' @param one_group_reps number of random starts for each numerical optimization in 1-component model
#' @param eps stopping tolerance for EM algoithm
#' @param max.it maximum number of iterations of EM algorithm
#' @param verbose print information as fitting progresses?
#'
#' @description This function generate a mixture of lognormal, gamma, and weibull distributions
#'
#' @details Data set x and number of components G are provided
#'
#' @return an ltmm model object
#'
#' @export
ltmm <- function(
    x, # data vector
    G, # number of components
    distributions, # density of each component (if length = 1, density of all components)
    trunc = NULL, # left truncation point (optional)
    EM_init_method = "emEM", # initialization method for EM algorithm
    EM_starts = 5, # number of random starts for initialization of EM algorithm. (only for G > 1)
    init_pars = NULL, # manually specified initial parameter values (for init_method=specified)
    init_pi = NULL, # manually specified initial component proportions (for init_method=specified)
    init_classes = NULL, # manually specified initial classes. will overwrite init_pars and init_pi
    one_group_reps = 50, # number of repetitions for each numerical optimization
    eps = 1e-6, # stopping tolerance for EM algoithm
    max.it = 1000, # maximum number of iterations for EM algorithm
    verbose = FALSE # print information as fitting progresses?
) {
    if (G != 1 && length(distributions) == 1)
        distributions = rep(distributions, G)
    if (!all(distributions %in% c("normal", "lognormal", "gamma", "weibull", "gompertz", "burr", "inv_burr")))
        stop("Invalid choice of density.")

    # create density functions
    f_list <- lapply(distributions, create_trunc_pdf, trunc=trunc)
    # record names
    attr(f_list, "names") <- distributions

    # package EM options
    EM_options <- list(eps = eps, max.it = max.it, verbose = verbose)

    if (G == 1) {

        if (is.null(init_pars))
            init_pars <- generate_initial_pars(x, distributions)[[1]]

        Pars <- single_group_mle(x, f_list[[1]], distributions, trunc, init_pars, one_group_reps, verbose)
        ll <- log_likelihood(x, f_list, Pars, Pi = 1)
        npars <- length(unlist(Pars))
        n <- length(x)

        # make single-component fit match multi-component fits:

        EM_result <- list(
            G = 1, Pi = 1, Pars = Pars,
            ll = ll,
            bic = -2 * ll + npars * log(n),
            aic = -2 * ll + npars * 2,
            id = rep(1, n),
            iter = 1,
            npars = npars,
            ll.history = ll
        )
    } else {

        if (EM_init_method == "spec")
            initial <- list(Pars = init_pars, Pi = init_pi)
        else
            initial <- initialize_EM(x, G, EM_starts, distributions, f_list, EM_init_method, init_classes, verbose)

        if (verbose)
            cat("Beginning EM iterations.\n\n")

        EM_result <- EM(x, initial$Pars, initial$Pi, f_list, EM_options)
    }

    if (verbose)
        cat("\nModel fitting completed.\n\n")

    # create fitted pdf
    fitted_pdf <- function(X) sapply(X, function(y)
        sum(sapply(1:G, function(g)
            EM_result$Pi[g]*f_list[[g]](y, EM_result$Pars[[g]]))))

    # create fitted cdf
    cdf_list <- lapply(distributions, create_trunc_cdf, trunc=trunc)

    fitted_cdf <- function(X) sapply(X, function(y) {
        sum(sapply(1:G, function(g)
            EM_result$Pi[g]*cdf_list[[g]](y, EM_result$Pars[[g]])))
    })

    VaR <- create_VaR_func(x, trunc, fitted_pdf, fitted_cdf)
    ES <- create_ES_func(x, trunc, fitted_pdf, fitted_cdf, VaR)

    structure(
        c(list(x=x,
               distributions=distributions,
               trunc=trunc,
               fitted_pdf = fitted_pdf,
               fitted_cdf = fitted_cdf,
               VaR = VaR,
               ES = ES),
          EM_result),
        class = "ltmm")
}

# Model Fitting ====

EM <- function(x, Pars, Pi, f_list, options) {
    # Main expectation-maximization algorithm loop.
    n <- length(x)
    G <- length(Pi)
    iter <- 0
    ll.history <- numeric(options$max.it + 2)

    ll.old <- log_likelihood(x, f_list, Pars, Pi)
    ll.history[1] <- ll.old

    if (options$verbose) print_status(iter, Pi, Pars, ll.old)

    repeat {
        z <- E_step(x, f_list, Pars, Pi)
        W <- M_step(x, f_list, z, Pars)

        Pi <- W$Pi
        Pars <- W$Pars

        ll <- log_likelihood(x, f_list, Pars, Pi)
        ll.history[iter+2] <- ll

        iter <- iter + 1
        if (options$verbose) print_status(iter, Pi, Pars, ll)

        if (abs((ll - ll.old) / ll) < options$eps) break
        if (iter >= options$max.it) break

        ll.old <- ll

    }

    id.est <- apply(z, 1, which.max)
    npars <- length(unlist(Pars)) + (G - 1)
    aic <- -2 * ll + npars * 2
    bic <- -2 * ll + npars * log(n)

    list(
        G = G,
        Pi = Pi,
        Pars = Pars,
        ll = ll,
        bic = bic,
        aic = aic,
        id = id.est,
        iter = iter,
        npars = npars,
        ll.history = ll.history[1:iter]
    )
}


E_step <- function(x, f_list, Pars, Pi) {
    # Expectation-step: update group membership weights (z)
    n <- length(x)
    G <- length(Pi)
    z <- matrix(NA, n, G)

    for (i in 1:n) {
        for (g in 1:G) {
            z[i,g] <- Pi[g] * f_list[[g]](x[i], Pars[[g]])
            if (is.nan(z[i,g])) z[i,g] <- 0
        }
        z[i,] <- z[i,] / sum(z[i,])
    }

    z
}


M_step <- function(x, f_list, z, Pars) {
    # M-step: Update parameter values conditional on z.
    n <- length(x)
    G <- ncol(z)
    Pi <- numeric(G)

    for (g in 1:G) {
        Pi[g] <- mean(z[,g])
        max_Q <- suppressWarnings(optim(Pars[[g]], fn=Q, f=f_list[[g]], x=x, z=z[,g]))
        Pars[[g]] <- max_Q$par
    }

    list(Pi = Pi, Pars = Pars)
}

generate_initial_pars <- function(x, distributions, partition = NULL) {
    # Generates initial parameter estimates for the lognormal, gamma and Weibull distributions,
    # for a set of observations x.
    # The lognormal starting values are based on the moments of a non-truncated lognormal
    # distribution and will thus be biased for a non-truncated distribution.
    # The values for the gamma and Weibull are based on the starting values for the non-truncated
    # MLE from the MASS::fitdistr function.

    G <- length(distributions)
    if (G == 1) {
        partition <- rep(1, length(x))
    }

    Pars <- list()

    for (g in 1:G) {
        density_name <- distributions[g]
        x_g <- x[partition == g]
        if (density_name == "lognormal") {
            # initial values are non-truncated MLE
            lx <- log(x_g)

            Pars[[g]] <- c(meanlog = mean(lx), sdlog = sd(lx))
        }

        else if (density_name == "gamma") {
            # initial paramater values analogous to MASS::fitdistr function for non-truncated MLE
            m <- mean(x_g)
            v <- var(x_g)

            Pars[[g]] <- c(shape = m^2/v, scale  = v/m)
        }
        else if (density_name == "weibull") {
            # initial paramater values same as MASS::fitdistr function for non-truncated MLE
            lx <- log(x_g)
            m <- mean(lx)
            v <- var(lx)
            shape <- 1.2/sqrt(v)
            scale <- exp(m + 0.572/shape)

            Pars[[g]] <- c(shape = shape, scale = scale)
        }
    }

    Pars
}

generate_random_initial_pars <- function(m, density_name, base) {
    # Generates a set of m random initial starting loosely based a starting set of values (base).
    pars <- matrix(NA, m, 2)
    if (density_name == "lognormal") {
        for (i in 1:m)
            pars[i,] <- c(log(rexp(1, 1/exp(base[1]))), rexp(1, 1/base[2]))
    } else if (density_name %in% c("gamma", "weibull"))
        for (i in 1:m)
            pars[i,] <- c(rexp(1, 1/base[1]), rexp(1, 1/base[2]))

        rbind(base, pars)
}


single_group_mle <- function(x, density_function, density_name, trunc, init_pars, one_group_reps, verbose = FALSE) {
    # Finds MLE of a single component LTMM. Repeatedly performs optimization with optim using perturbed
    # starting parameters, and selects the best overall optimization.

    if (one_group_reps > 0) {
        random_initial_pars <- generate_random_initial_pars(one_group_reps, density_name, init_pars)
        attempts <- suppressWarnings(
            t(apply(random_initial_pars, 1,
                    function(pars) {
                        tryCatch(
                            unlist(optim(pars, fn = Q, f = density_function, x = x)[c('par', 'value')]),
                            error = function(e)
                                c(NA, NA, NA))
                    }
            ))
        )
        successful_attempts <- attempts[!is.na(attempts[,3]), ]

        best.ind <- which.min(successful_attempts[,3])
        pars <- successful_attempts[best.ind, 1:2]

        attemptsSummary <- cbind(random_initial_pars, attempts)
        colnames(attemptsSummary) = c("par1_i", "par2_i", "par1_f", "par2_f", "Q")

    } else
        pars <- suppressWarnings(optim(init_pars, fn=Q, f=density_function, x=x)$par)

    if (verbose) {
        cat(one_group_reps+1, "optimization attempts.", nrow(successful_attempts), "successful.\n\n")
        print(attemptsSummary)
    }

    list(pars)
}

# Function for optimization: negative conditional expectation of the complete log likelihood
Q <- function(pars, f, x, z = NULL) {
    if (is.null(z))
        -sum(f(x, pars, as.log = TRUE))
    else
        -sum(z * f(x, pars, as.log = TRUE))
}

log_likelihood <- function(x, f_list = NULL, Pars, Pi, trunc, distributions) {
    G <- length(Pi) # number of groups
    n <- length(x) # number of observations

    if (is.null(f_list))
        f_list <- lapply(distributions, create_trunc_pdf, trunc = trunc)

    sum(sapply(1:n, function(i)
        log(sum(sapply(1:G, function(g)
            Pi[g] * f_list[[g]](x[i], Pars[[g]])))))
    )
}

generate_partition <- function(x, G) {
    # generates G random partitions of x based on random centers
    n <- length(x)
    centers <- sample(x, G, replace=FALSE)
    dis <- matrix(0, n, G)
    for(g in 1:G) dis[,g] <- abs(x - centers[g])

    apply(dis, 1, which.min)
}

initialize_EM <- function(x, G, EM_starts, distributions, f_list, method, init_classes, verbose){
    # creates random partitions and corresponding parameter values to initialize the EM algorithm
    # only used when G > 1

    attempts <- list()

    if (verbose)
        cat(paste0("Initializing using ", method, " initialization.\n"))

    for (i in 1:EM_starts) {
        # populate EM_starts attempts
        if (verbose)
            cat(paste("\nInitialization attempt", i, "\n\n"))
        if (is.null(init_classes))
            random_partition <- generate_partition(x, G)
        else
            random_partition <- init_classes

        init_pi <- table(random_partition) / length(random_partition)
        init_pars <- generate_initial_pars(x, distributions, random_partition)

        EM_options <- list(eps = 1e-1, max.it = 100, verbose = verbose)

        if (method == "RndEM") {
            attempts[[i]] <- list(Pars = init_pars, Pi = init_pi, ll = log_likelihood(x, f_list, init_pars, init_pi))
        } else if (method == "emEM") {
            attempts[[i]] <- suppressWarnings(
                tryCatch(
                    expr = EM(x, init_pars, init_pi, f_list, EM_options)[c("Pars", "Pi", "ll")],
                    error = function(e) list(ll = -Inf)
                )
            )
        } else
            stop("Invalid initialization method.")
    }

    log_likelihoods <- sapply(attempts, getElement, "ll")
    if (all(is.infinite(log_likelihoods)))
        stop("No initializations were successful for the EM algorithm.")

    best_attempt <- attempts[[which.max(log_likelihoods)]]

    if (verbose)
        cat("\nIntialization successful.\n\n")

    best_attempt
}


# Take density name & produce pdf as function object
select_pdf <- function(density_name) {
    switch(density_name,
           "normal"    = function(x, pars, as.log=FALSE) dnorm(x, mean = pars[1], sd = sqrt(pars[2]), log=as.log),
           "lognormal" = function(x, pars, as.log=FALSE) dlnorm(x, meanlog = pars[1], sdlog = pars[2], log=as.log),
           "gamma"     = function(x, pars, as.log=FALSE) if(!all(pars > 0)) return(Inf)
           else dgamma(x, shape = pars[1], scale = pars[2], log=as.log),
           "weibull"   = function(x, pars, as.log=FALSE) {
               # dweibull sometimes returns NaN instead of zero --- annoying
               y <- suppressWarnings(dweibull(x, shape = pars[1], scale = pars[2], log=as.log))
               if (is.nan(y)) if(as.log) -Inf else 0 else y
           },
           "gompertz"  = function(x, pars, as.log=FALSE) flexsurv::dgompertz(x, shape = pars[1], rate = pars[2], log=as.log),
           "burr"      = function(x, pars, as.log=FALSE) actuar::dburr(x, shape1 = pars[1], shape2 = pars[2], scale = pars[3], log=as.log),
           "inv_burr"  = function(x, pars, as.log=FALSE) actuar::dinvburr(x, shape1 = pars[1], shape2 = pars[2], scale = pars[3], log=as.log)
    )
}

# Take density name & produce cdf as function object
select_cdf <- function(density_name) {
    switch(density_name,
           "normal"    = function(x, pars) pnorm(x, mean = pars[1], sd = sqrt(pars[2])),
           "lognormal" = function(x, pars) plnorm(x, meanlog = pars[1], sdlog = pars[2]),
           "gamma"     = function(x, pars) pgamma(x, shape = pars[1], scale = pars[2]),
           "weibull"   = function(x, pars) pweibull(x, shape = pars[1], scale = pars[2]),
           "gompertz"  = function(x, pars) flexsurv::pgompertz(x, shape = pars[1], rate = pars[2]),
           "burr"      = function(x, pars) actuar::pburr(x, shape1 = pars[1], shape2 = pars[2], scale = pars[3]),
           "inv_burr"  = function(x, pars) actuar::pinvburr(x, shape1 = pars[1], shape2 = pars[2], scale = pars[3])
    )
}

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

create_ES_func <- function(x, trunc, fitted_pdf, fitted_cdf, VaR) {
    function(p = 0.05) {
        if (p < 0 | p > 1)
            stop("p must be between 0 and 1.")
        q0 <- VaR(p)
        fx <- function (x) fitted_pdf(x)*x

        pracma::quadinf(fx, q0, Inf)$Q / pracma::quadinf(fitted_pdf, q0, Inf)$Q
    }
}


# Display & Analysis ====

# Print update for each iteration
print_status <- function(iter, Pi, Pars, ll) {
    cat("\tIteration", iter, "\n")
    cat("\t\tpi =", Pi, "\n")
    cat("\t\tPars =", "\n")
    sapply(Pars,function(x) cat("\t\t", x,'\n'))
    cat("\t\tll =", ll, "\n")
}

# Plot method for ltmm objects.
# Displays histogram & fitted pdf.
plot.ltmm <- function(fit, bins = 40, xlim = NULL, ylim = NULL) {
    with(fit,
         {

             if (!is.null(xlim))
                 x = x[x > xlim[1] & x < xlim[2]]
             f_list <- lapply(distributions, create_trunc_pdf, trunc=trunc)
             fitted_dens <- function(y)
                 sum(sapply(1:G, function(g) Pi[g]*f_list[[g]](y, Pars[[g]])))

             hist(x, bins, probability=TRUE, main="Fitted Density", ylim = ylim)
             t <- seq(min(x),max(x),length.out=200)
             lines(t, sapply(t, fitted_dens), lwd=2)
             for (g in 1:G) lines(t, sapply(t, function(y) Pi[g]*f_list[[g]](y, Pars[[g]])),col=g+1)

             legend_names = sapply(1:G, function(g) paste("Group",g))
             legend("topright", legend = legend_names, fill=(1:G)+1, border=F)
         }
    )
}

summary.ltmm <- function(fit) {
    cat("Model:", fit$distributions, "| BIC =", fit$bic, "| AIC =", fit$aic, "| log-likelihood =", fit$ll)
}

summary.ltmmCombo <- function(w) {
    cat("====", w$G, "COMPONENT COMBINATIONS OF", toupper(w$distributions), "DISTRIBUTIONS ====\n\n")
    cat("Best model by BIC:\n")
    with(w$best.bic.fit, cat("Model:", distributions, "| BIC =", bic, "| AIC =", aic, "| log-likelihood =", ll, "\n\n"))
    cat("Best model by AIC:\n")
    with(w$best.aic.fit, cat("Model:", distributions, "| BIC =", bic, "| AIC =", aic, "| log-likelihood =", ll, "\n\n"))
    print(w$summary_table)

    invisible(NULL)
}

createLtmmObj <- function(x, distributions, trunc, Pars, Pi, npars = NULL) {

    # creates an ltmm object from parameters. Useful for comparing models produced using ltmm
    # to models found in other ways

    # can use npars to overwrite the number of free parameters, if the model has additional
    # constraints

    G <- length(Pi)
    n <- length(x)
    if (G != 1 && length(distributions) == 1)
        distributions = rep(distributions, G)

    f_list <- lapply(distributions, create_trunc_pdf, trunc=trunc)

    ll <- log_likelihood(x = x, Pars = Pars, Pi = Pi, trunc = trunc, distributions = distributions)

    # create fitted pdf
    fitted_pdf <- function(X) sapply(X, function(y)
        sum(sapply(1:G, function(g)
            Pi[g]*f_list[[g]](y, Pars[[g]]))))

    # create fitted cdf
    cdf_list <- lapply(distributions, create_trunc_cdf, trunc=trunc)

    fitted_cdf <- function(X) sapply(X, function(y) {
        sum(sapply(1:G, function(g)
            Pi[g]*cdf_list[[g]](y, Pars[[g]])))
    })

    VaR <- create_VaR_func(x, trunc, fitted_pdf, fitted_cdf)
    ES <- create_ES_func(x, trunc, fitted_pdf, fitted_cdf, VaR)

    if (is.null(npars))
        npars <- length(unlist(Pars))
    aic <- -2 * ll + npars * 2
    bic <- -2 * ll + npars * log(n)

    id <- sapply(x, function(x0) {
        which.max(sapply(1:G, function(i) f_list[[i]](x0, Pars[[i]])))
    })

    structure(
        list(x=x,
             distributions = distributions,
             trunc = trunc,
             fitted_pdf = fitted_pdf,
             fitted_cdf = fitted_cdf,
             VaR = VaR,
             ES = ES,
             G = G,
             Pi = Pi,
             Pars = Pars,
             ll = ll,
             bic = bic,
             aic = aic,
             id = id,
             iter = 0,
             npars = npars,
             ll.history = ll ),
        class = "ltmm")

}

# Simulation Tools ====

selectGenerationFunction <- function(distribution, pars) {
    switch(distribution,
           "lognormal" = function(x) rlnorm(1, meanlog = pars[1], sdlog = pars[2]),
           "gamma"     = function(x) rgamma(1, shape = pars[1], scale = pars[2]),
           "weibull"   = function(x) rweibull(1, shape = pars[1], scale = pars[2])
    )
}

createSamplingFunction <- function(distributions, trunc, Pars, Pi) {

    G <- length(Pi)
    if (G != 1 && length(distributions) == 1)
        distributions = rep(distributions, G)

    rdist <- mapply(selectGenerationFunction, distributions, Pars)

    function(n) {
        output <- numeric(n)

        for (i in 1:n) {
            if (G == 1)
                g <- 1
            else
                g <- min(which(runif(1) < cumsum(Pi)))

            x <- -Inf
            while(x < trunc)
                x <- rdist[[g]]()

            output[i] <- x
        }

        output
    }

}













