## #########################################
## Simulation 1: Foerster 2016 like setup ##
## #########################################
## Simulate different number of interactions (balanced) and compare to ML approach.
## Perform prior sensitivity check by refitting using different standard deviations 
## for the sigma prior.

## Colors:
cols <- c(rgb(0.1, 0.5, 0.1, alpha = 0.6), rgb(0.2, 0.2, 0.2, alpha = 0.6))

## Estimation:
set.seed(123)
starting_scores <- rnorm(n = 44, mean = 0, sd = 2.62)
n_interaction <- c(1500, 1000, 500, 250)
## Comment out (until including 'save(...)'-commands to save time (after first time having saved the results)):
ml_starting_scores <- bayes_starting_scores <-
  bayes_starting_scores_sd2 <- vector("list", length(n_interaction))
for (i in 1:length(n_interaction)) {
  set.seed(1604)
  X_sim <- simulate_foerster_data(starting_scores = starting_scores,
                                  n_interaction = n_interaction[i])
  Ai <- apply(X_sim, MAR = 1, FUN = function(x){which(x == 1)})
  Bi <- apply(X_sim, MAR = 1, FUN = function(x){which(x == -1)})
  fit_dat <- list(N = nrow(X_sim), K = ncol(X_sim), Ai = Ai, Bi = Bi,
                  y = rep(1, nrow(X_sim)), diff_f = 1, presence = 1 + 0 * X_sim)
  chains <- 1; thin <- 1; iter <- 1500; warmup <- 500
  iter_used <- iter - warmup
  m <- stan(file = 'elo_score_USE_THIS.stan', data = fit_dat,
            iter = iter*thin, chains = chains, thin = thin, warmup = warmup*thin,
            control = list(adapt_delta = 0.95))
  bayes_starting_scores[[i]] <- colMeans(extract(m)[["EloStart"]])
  m <- stan(file = 'elo_score_USE_THIS_sd2.stan', data = fit_dat,
            iter = iter*thin, chains = chains, thin = thin, warmup = warmup*thin,
            control = list(adapt_delta = 0.95))
  bayes_starting_scores_sd2[[i]] <- colMeans(extract(m)[["EloStart"]])
  foerster_result_model3 <- optim(par = c(rep(0, ncol(X_sim)), log(1)),
                                  fn = logLik, X = X_sim, show_k = TRUE, method = 'BFGS',
                                  control = list(maxit = 10000, reltol = 1e-10), hessian = TRUE)
  ml_starting_scores[[i]] <- foerster_result_model3$par[-length(foerster_result_model3$par)]
}
save(ml_starting_scores, file = paste0("sim1_ml_starting_scores_", today, ".RData"))
save(bayes_starting_scores, file = paste0("sim1_bayes_starting_scores_", today, ".RData"))
save(bayes_starting_scores_sd2, file = paste0("sim1_bayes_starting_scores_", today, ".RData"))
## This will save time:
# load(file = "sim1_ml_starting_scores_", today, ".RData"))
# load(file = "sim1_bayes_starting_scores_", today, ".RData"))
# load(file = "sim1_bayes_starting_scores_sd2_", today, ".RData"))

## Plot results:
pdf(paste0("simulation1_", today, ".pdf"), height = 5, width = 8)
par(mfrow = c(2, 4), mar = c(4.1, 4.1, 3, 2))
for (i in 1:length(n_interaction)) {
  plotting_range <- range(c(bayes_starting_scores[[i]], ml_starting_scores[[i]], 
                            starting_scores))
  plot(bayes_starting_scores[[i]], starting_scores, pch = 16, bty = "n", 
       yaxt = "n", xlab = "Estimated start rating", ylab = "True start rating", 
       xlim = plotting_range, ylim = plotting_range, 
       main = paste0("No. of interactions: ", n_interaction[i], "\nFull view"), 
       col = cols[1], cex.main = 0.9)
  axis(2, las = 2); abline(a = 0, b = 1, lty = 2)
  points(ml_starting_scores[[i]], starting_scores, col = cols[2], pch = 16)
  legend("bottomright", pch  = 16, col = rev(cols), legend = c("ML", "Bayes"), 
         bty = "n", cex = 0.6, pt.cex = 1)
}
for (i in 1:length(n_interaction)) {
  plotting_range <- range(-8, 8)
  plot(bayes_starting_scores[[i]], starting_scores, pch = 16, bty = "n", yaxt = "n", 
       xlab = "Estimated start rating", ylab = "True start rating", 
       xlim = plotting_range, ylim = plotting_range, 
       main = paste0("No. of interactions: ", n_interaction[i], "\nEnlarged view of [-8, 8]"), 
       col = cols[1], cex.main = 0.9)
  axis(2, las = 2); abline(a = 0, b = 1, lty = 2)
  points(ml_starting_scores[[i]], starting_scores, col = cols[2], pch = 16)
  ## Those arrows were too tiny for publication figure (but maybe still helpful on a screen?):
  # arrows(x0 = bayes_starting_scores[[i]], y0 = starting_scores, 
  #        x1 = bayes_starting_scores_sd2[[i]], y1 = starting_scores, length = 0.015, lwd = 0.5)
  for (j in 1:length(bayes_starting_scores[[i]])) {
    lines(c(bayes_starting_scores[[i]][j], bayes_starting_scores_sd2[[i]][j]), 
          c(starting_scores[j], starting_scores[j]), col = cols[1])
  }
  points(bayes_starting_scores_sd2[[i]], starting_scores, col = cols[1], pch = 4)
  legend("topleft", pch  = c(16, 16, 4), col = cols[c(2, 1, 1)], cex = 0.6, pt.cex = 1,
         legend = c("ML", "BI (smaller prior variance)", "BI (larger prior variance)"), bty = "n")
}
dev.off()
shell.exec(paste0("simulation1_", today, ".pdf"))