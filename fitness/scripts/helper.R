test_for_overdispersion <- function(model){ 
  chisq <- sum(resid(model,type='pearson')^2)
  df_resid <- df.residual(model)
  residual_deviance = chisq/df_resid
  chi2test <- pchisq(chisq, df_resid, lower.tail=FALSE)
  return(c("residual_deviance" = residual_deviance, "chi2test" = chi2test))
}



prediction <- function(model, pd, log = FALSE) {
  # Generate predicted values on the response scale
  pd$NumFliesPerPlant <- predict(model, pd, type = "response")
  
  # Generate predicted values on the link scale with standard errors
  pred <- predict(model, pd, type = "link", se.fit = TRUE)
  
  # Extract the fitted values and standard errors
  fit_link <- pred$fit
  se_link <- pred$se.fit
  
  # Calculate the 95% CI on the link scale
  upper_ci_link <- fit_link + 1.96 * se_link
  lower_ci_link <- fit_link - 1.96 * se_link
  
  if(log) {
    fit <- exp(fit_link)
    upper_ci <- exp(upper_ci_link)
    lower_ci <- exp(lower_ci_link)
  } else {
    fit <- fit_link
    upper_ci <- upper_ci_link
    lower_ci <- lower_ci_link
  }
  
  # Combine into a data frame
  results <- data.frame(
    pd,
    fit = fit,
    upper_ci = upper_ci,
    lower_ci = lower_ci
  )
  
  # Return the results
  return(results)
}




plot_prediction <- function(predictons, ymax){ 
  ggplot(predictons, aes(x = PhenotypeCageSpecies, y = fit, color = PopSpecies, group = PopSpecies)) +
    geom_point(size = 3, stroke = 2) +
    ylim(0, ymax) + 
    geom_line() + 
    geom_errorbar(aes(ymin = lower, ymax = upper, color = PopSpecies), width = 0.2, size = 1) +
    labs(title = "Predictions with Confidence Intervals",
         x = "Phenotype Cage Species",
         y = "Number of Flies Per Plant") + 
    facet_wrap(~PopSpecies)
  }


sim_neg_binomial <- function(model, n, covtable){
  fitted_vals <- round(predict(model, covtable, type = "response"))
  theta <- model$theta
  nsim = n
  msim <- matrix(ncol = dim(covtable)[1], nrow = nsim)
  for(i in 1:nsim){
    covtable$sim <- rnegbin(n = nrow(covtable), mu = fitted_vals, theta = theta)
    smod <- glm.nb(sim ~ trays + PopSpecies + PhenotypeCageSpecies + PopSpecies:PhenotypeCageSpecies, data = covtable)
    test_for_overdispersion(smod)
    msim[i, ] <- predict(smod, covtable, type = "response")
  }
  return(msim)
}


# From https://stackoverflow.com/questions/72820236/comparing-all-factor-levels-to-the-grand-mean-can-i-tweak-contrasts-in-linear-m
# https://stackoverflow.com/questions/41032858/lm-summary-not-display-all-factor-levels
ContrSumMat <- function (fctr, sparse = FALSE) {
  if (!is.factor(fctr)) stop("'fctr' is not a factor variable!")
  N <- nlevels(fctr)
  Cmat <- contr.sum(N, sparse = sparse)
  dimnames(Cmat) <- list(levels(fctr), seq_len(N - 1))
  Cmat
}
ContrSumMat(as.factor(tab$PhenotypeCageSpecies))


edit_plot_data <- function(x){
  plotd <- x
  plotd$PopSpecies <- sub(x = plotd$PopSpecies, "^B$", "Barbarea")
  plotd$PopSpecies <- sub(x = plotd$PopSpecies, "^T$", "Turritus")
  plotd$PopSpecies <- sub(x = plotd$PopSpecies, "^M$", "Mixture")
  plotd$PhenotypeCageSpecies <- sub(x = plotd$PhenotypeCageSpecies, "^B$", "Tested\n on Barbarea")
  plotd$PhenotypeCageSpecies <- sub(x = plotd$PhenotypeCageSpecies, "^T$", "Tested\n on Turritus")
  plotd$PhenotypeCageSpecies <- sub(x = plotd$PhenotypeCageSpecies, "^M$", "Tested\n on Mixture")
  plotd$lines <- paste(plotd$PopNumber, plotd$PopSpecies)
  plotd$PopNumber <- sub(x = plotd$PopNumber, "^3$", "Source Population 1")
  plotd$PopNumber <- sub(x = plotd$PopNumber, "^4$", "Source Population 2")
  return(plotd)
}

glm_sim_nb_slim <- function(x, n_sims, dframe, max_iterations = 10000){
  dat <- dframe
  sim_results <- vector("list", n_sims)
  valid_sims <- 0
  iteration <- 1
  while (valid_sims < n_sims) {
    tryCatch({
      withCallingHandlers({
        form <- paste(deparse(eval(x$call[[2]])), collapse = " ")
        form2 <- sub(x = form, pattern = "numFlies", replacement = "y_sim")
        mu <- predict(x, newdata = dat, type = "response")
        theta <- x$theta
        dat$y_sim <- rnegbin(n = length(mu), mu = mu, theta = theta)
        mod_sim <- glm.nb(y_sim ~ 
                            PopSpecies + 
                            PhenotypeCageSpecies + 
                            PopSpecies:PhenotypeCageSpecies + 
                            offset(log(trays)),
                          contrasts = list(PopSpecies = contrast_matrix, PhenotypeCageSpecies = contrast_matrix),
                          data = dat)
      }, warning = function(w) {
        stop("Warning occurred")
      })
      sim_results[[valid_sims + 1]] <- predict(mod_sim, dat, type = "response")
      valid_sims <- valid_sims + 1
    }, error = function(e) {
    })
    iteration <- iteration + 1
    if (iteration > max_iterations) {
      warning(paste("Reached maximum iterations. Only", valid_sims, "valid simulations were completed."))
      break
    }
  }
  sim_results <- sim_results[!sapply(sim_results, is.null)]
  sim_results_matrix <- do.call(cbind, sim_results)
  dat$lower_ci <- apply(sim_results_matrix, 1, quantile, probs = 0.025)
  dat$upper_ci <- apply(sim_results_matrix, 1, quantile, probs = 0.975)
  print(sim_results_matrix)
  return(dat) 
}


