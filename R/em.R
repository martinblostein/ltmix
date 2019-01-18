# These functions are concerned with fitting an individual model via the EM algorithm

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

# Print update for each iteration
print_status <- function(iter, Pi, Pars, ll) {
  cat("\tIteration", iter, "\n")
  cat("\t\tpi =", Pi, "\n")
  cat("\t\tPars =", "\n")
  sapply(Pars,function(x) cat("\t\t", x,'\n'))
  cat("\t\tll =", ll, "\n")
}
