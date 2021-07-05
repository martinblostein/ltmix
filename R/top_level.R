# These functions are the top-level interface functions that are exposed to the user of the package.

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
#' @param save_each_fit save each model as it is produced, in a time-stamped directory (safer)
#'
#' @description This function fits a family of finite mixture models using every combination of the left-truncated
#'              lognormal, gamma, and weibull distributions.
#'
#' @references Blostein, Martin & Miljkovic, Tatjana. (2019). On modeling left-truncated loss data using mixtures of distributions. Insurance Mathematics and Economics. 85. 35-46. 10.1016/j.insmatheco.2018.12.001.
#'
#' @examples
#' \donttest{
#' x <- secura$Loss
#'
#' fits_GL <- ltmmCombo(x, G = 2, distributions = c('gamma', 'lognormal'), trunc = 1.2e6)
#' summary(fits_GL)
#' }
#'
#' @return An ltmmCombo model object, with the following properties:
#' \describe{
#'   \item{x}{Copy of the input data}
#'   \item{distributions}{The selected distributions}
#'   \item{combos}{List of all combinations of distributions considered}
#'   \item{all.fits}{List of all ltmm fit objects}
#'   \item{all.bic}{Vector of BIC values for each model}
#'   \item{best.bic.fit}{The best ltmm fit by BIC}
#'   \item{best.bic}{The best BIC value of all fits}
#'   \item{best.bic.combo}{The combination of distributions used for the best fit by BIC}
#'   \item{all.aic}{Vector of AIC value for each model}
#'   \item{best.aic.fit}{The best ltmm fit by AIC}
#'   \item{best.aic}{The best AIC value of all fits}
#'   \item{best.aic.combo}{The combination of distributions used for the best fit by AIC}
#'   \item{all.ll}{Vector of log-likelihood value for each model}
#'   \item{summary_table}{Table summarizing the AIC, BIC, LL, and risk measures for each fitted model}
#' }
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
#' @description This function generates a mixture model combining left-truncated lognormal,
#'              gamma, and weibull distributions
#' @examples
#' \donttest{
#' x <- secura$Loss
#'
#' fit <- ltmm(x, G = 2, distributions = c('gamma', 'gamma', 'weibull'), trunc = 1.2e6)
#'
#' summary(fit)
#' plot(fit)
#' }
#'
#' @return An ltmm model object, with the following properties:
#' \describe{
#'   \item{x}{Copy of the input data}
#'   \item{distributions}{The selected distributions}
#'   \item{trunc}{The left truncation value, if specified}
#'   \item{fitted_pdf}{The probability density function of the fitted model}
#'   \item{fitted_cfd}{The cumulative density function of the fitted model}
#'   \item{VaR}{The value-at-risk of the fitted model (function with p taken as onl yargument)}
#'   \item{ES}{The expected shortfall of the fitted model (function with p taken as onl yargument)}
#'   \item{G}{The number of components in the model}
#'   \item{Pi}{The estimated probabilites of component membership}
#'   \item{Pars}{The estimated model parameters}
#'   \item{ll}{The log-likelihood of the fitted model}
#'   \item{bic}{The BIC of the fitted model}
#'   \item{aic}{The AIC of the fitted model}
#'   \item{id}{The MAP component membership for each observation}
#'   \item{iter}{The number of iterations until convergence for the EM algorithm}
#'   \item{npars}{The total number of model parameters for the fitted model}
#'   \item{ll.history}{The value of log-likelihood at each iteration of the EM algorithm}
#' }
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

#' Create an ltmm model object given data and parameters
#'
#' @param x data vector
#' @param distributions densities to combine
#' @param trunc left truncation point (optional)
#' @param Pars list of length G of parameter values
#' @param Pi vector of length G of component proportions
#' @param npars Can optionally be used to overwrite the number of free parameters (used in the
#'        calculation of AIC & BIC), if the model has additional constraints
#'
#' @description This function is useful for omparing models produced using the ltmix package
#'               to models fit using other, or for computing fit criteria and risk measures for
#'               a known set of parameters.
#'
#' @return An ltmm model object
#'
#' @export
createLtmmObj <- function(x, distributions, trunc, Pars, Pi, npars = NULL) {

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
