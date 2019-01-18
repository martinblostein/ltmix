# These functions define S3 class methods for the classes outputted by the user-facing functions,
# ltmm and ltmmCombo.

# Plot method for ltmm objects: displays histogram & fitted pdf
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

# Summary method for ltmm objects: displays the model fitting criteria for the given model
summary.ltmm <- function(fit) {
  cat("Model:", fit$distributions, "| BIC =", fit$bic, "| AIC =", fit$aic, "| log-likelihood =", fit$ll)
}

# Summary method for ltmmCombo objects: displays a list of the best models by each model fitting
# criterion and a summary table
summary.ltmmCombo <- function(w) {
  cat("====", w$G, "COMPONENT COMBINATIONS OF", toupper(w$distributions), "DISTRIBUTIONS ====\n\n")
  cat("Best model by BIC:\n")
  with(w$best.bic.fit, cat("Model:", distributions, "| BIC =", bic, "| AIC =", aic, "| log-likelihood =", ll, "\n\n"))
  cat("Best model by AIC:\n")
  with(w$best.aic.fit, cat("Model:", distributions, "| BIC =", bic, "| AIC =", aic, "| log-likelihood =", ll, "\n\n"))
  print(w$summary_table)

  invisible(NULL)
}
