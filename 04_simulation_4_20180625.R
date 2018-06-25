## ############################################
## Mortality-Instability-Recovery Simulation ##
## ############################################
## Simulate 3000 interactions (balanced) with 
## an unstable period of 1000 interactions in the middle

cols <- viridis(n = 22, option = "A", alpha = 0.5, begin = 0, end = 0.8)
cols_without_alpha <- viridis(n = 22, option = "A", begin = 0, end = 0.8)

set.seed(1604)
starting_scores1 <- seq(-6, 6, length = 22) 
starting_scores2a <- sample(starting_scores1[-c(21, 22)])  
starting_scores2b <- sample(starting_scores1[-c(21, 22)])  
starting_scores2c <- sample(starting_scores1[-c(21, 22)])  
starting_scores2d <- sample(starting_scores1[-c(21, 22)])  
starting_scores2e <- sample(starting_scores1[-c(21, 22)])  
starting_scores2f <- sample(starting_scores1[-c(21, 22)])  
starting_scores2g <- sample(starting_scores1[-c(21, 22)])  
starting_scores2h <- sample(starting_scores1[-c(21, 22)])  
starting_scores2i <- sample(starting_scores1[-c(21, 22)])  
starting_scores2j <- sample(starting_scores1[-c(21, 22)])  
starting_scores3 <- starting_scores1[-c(21, 22)]  

chains <- 1; thin <- 10; iter <- 1500; warmup <- 500
iter_used <- iter - warmup

repeats <- 2
set.seed(1604)
fit_list <- X_list <- vector("list", repeats)
for (rep_index in 1:repeats) {
  X_sim1 <- simulate_mortality_instability_recovery_data(starting_scores = starting_scores1, 
                                                         n_interaction = 1000)
  X_sim2a <- simulate_mortality_instability_recovery_data(starting_scores = starting_scores2a, 
                                                          n_interaction = 100)
  X_sim2b <- simulate_mortality_instability_recovery_data(starting_scores = starting_scores2b, 
                                                          n_interaction = 100)
  X_sim2c <- simulate_mortality_instability_recovery_data(starting_scores = starting_scores2c, 
                                                          n_interaction = 100)
  X_sim2d <- simulate_mortality_instability_recovery_data(starting_scores = starting_scores2d, 
                                                          n_interaction = 100)
  X_sim2e <- simulate_mortality_instability_recovery_data(starting_scores = starting_scores2e, 
                                                          n_interaction = 100)
  X_sim2f <- simulate_mortality_instability_recovery_data(starting_scores = starting_scores2f, 
                                                          n_interaction = 100)
  X_sim2g <- simulate_mortality_instability_recovery_data(starting_scores = starting_scores2g, 
                                                          n_interaction = 100)
  X_sim2h <- simulate_mortality_instability_recovery_data(starting_scores = starting_scores2h, 
                                                          n_interaction = 100)
  X_sim2i <- simulate_mortality_instability_recovery_data(starting_scores = starting_scores2i, 
                                                          n_interaction = 100)
  X_sim2j <- simulate_mortality_instability_recovery_data(starting_scores = starting_scores2j, 
                                                          n_interaction = 100)
  X_sim3 <- simulate_mortality_instability_recovery_data(starting_scores = starting_scores3, 
                                                         n_interaction = 2000)
  X_sim <- rbind(X_sim1, 
                 cbind(X_sim2a, 0, 0), cbind(X_sim2b, 0, 0), cbind(X_sim2c, 0, 0), 
                 cbind(X_sim2d, 0, 0), cbind(X_sim2e, 0, 0), cbind(X_sim2f, 0, 0), 
                 cbind(X_sim2g, 0, 0), cbind(X_sim2h, 0, 0), cbind(X_sim2i, 0, 0), 
                 cbind(X_sim2j, 0, 0), 
                 cbind(X_sim3, 0, 0))
  X_list[[rep_index]] <- X_sim
  fit_dat <- list(N = nrow(X_sim), K = ncol(X_sim), X = X_sim)
  
  Ai <- apply(X_sim, MAR = 1, FUN = function(x){which(x == 1)})
  Bi <- apply(X_sim, MAR = 1, FUN = function(x){which(x == -1)})
  fit_dat <- list(N = nrow(X_sim), K = ncol(X_sim), Ai = Ai, Bi = Bi, 
                  y = rep(1, nrow(X_sim)), diff_f = 1)
  fit_dat$presence <- 1 + 0 * X_sim
  fit_dat$presence[1001:nrow(X_sim), ncol(X_sim) + c(-1, 0)] <- 0
  ## Comment out following command to save time (after first time having saved the results):
  fit_list[[rep_index]] <- stan(file = 'elo_score_USE_THIS.stan', data = fit_dat, 
                                iter = iter*thin, chains = chains, thin = thin, 
                                warmup = warmup*thin, control = list(adapt_delta = 0.95))
}
save(fit_list, file = paste0("sim4_fit_list_", today, ".RData"))
## This will save time:
# load(file = paste0("sim4_fit_list_", today, ".RData"))
## Plot results:
Elo_pA <- function(EloStart_logk, X, show_k = FALSE, presence){
  EloStart <- EloStart_logk[-length(EloStart_logk)]
  ## EloStart <- EloStart - mean(EloStart) ## added January 19, 2017
  k <- exp(EloStart_logk[length(EloStart_logk)])
  if (show_k) {
    cat(paste0(round(k, 3), paste(rep(" ", 20), collapse = " ")))
    cat("\r")
  }
  EloNow <- EloStart
  Elo <- matrix(nrow = nrow(X), ncol = length(EloStart), 0)
  colnames(Elo) <- colnames(X)
  pA <- rep(0, nrow(X))
  for (i in 1:nrow(X)) {
    A <- which(X[i, ] == 1)
    B <- which(X[i, ] == -1)
    EloNow <- EloNow - mean(EloNow[which(presence[i, ] == 1)])
    pA[i] <- pAfun(EloA = EloNow[A], EloB = EloNow[B])
    toAdd <- (1 - pA[i]) * k
    EloNow[A] <- EloNow[A] + toAdd
    EloNow[B] <- EloNow[B] - toAdd
    Elo[i, ] <- EloNow
  }
  return(list(pA = pA, Elo = Elo))
}

pdf(paste0("simulation4_", today, ".pdf"), height = 5, width = 7)
layout(mat = matrix(nrow = 2, ncol = 2, c(1, 3, 2, 4), byrow = T), 
       heights = c(0.65, 0.35), width = c(0.5, 0.5))
par(mar = c(4, 4, 1, 1))
for (rep_index in 1:repeats) {
  bayes_model3 <- fit_list[[rep_index]]
  X_sim <- X_list[[rep_index]]
  bayes_starting_scores <- summary(bayes_model3)$summary[1:ncol(X_sim)]
  bayes_k <- summary(bayes_model3)$summary[ncol(X_sim) + 1]
  elo_sim <- Elo_pA(EloStart_logk = c(bayes_starting_scores, log(bayes_k)), 
                    X = X_sim, presence = fit_dat$presence)$Elo
  plot(elo_sim[, 1], type = "n", ylim = range(c(elo_sim, starting_scores1, starting_scores3)), bty = "n", yaxt = "n", 
       xlab = "Interaction index", ylab = "Elo-rating", xlim = c(-100, 4100))
  bayes_k <- summary(bayes_model3)$summary[ncol(X_sim) + 1]
  main_title <- paste("k: ", round(bayes_k, 2), ", with 95% CI: [", 
                      paste0(unname(round(summary(bayes_model3)$summary[ncol(X_sim) + 1, c(4, 8)], 2)), 
                             collapse = ", "), "]", sep = "")
  title(main = main_title)
  axis(2, lab = FALSE, at = seq(-5, 5, by = 5))
  mtext(side = 2, at = seq(-5, 5, by = 5), text = seq(-5, 5, by = 5), cex = 1, line = 0.7, las = 2)
  lty_period1 <- rep(1, 22)
  lty_period2 <- c(rep(1, 20), 0, 0)
  points(rep(1, 22), starting_scores1, col = cols, pch = 16, cex = 0.95)
  points(rep(1050, 20), starting_scores2a, col = cols, pch = 16, cex = 0.95)
  points(rep(1150, 20), starting_scores2b, col = cols, pch = 16, cex = 0.95)
  points(rep(1250, 20), starting_scores2c, col = cols, pch = 16, cex = 0.95)
  points(rep(1350, 20), starting_scores2d, col = cols, pch = 16, cex = 0.95)
  points(rep(1450, 20), starting_scores2e, col = cols, pch = 16, cex = 0.95)
  points(rep(1550, 20), starting_scores2f, col = cols, pch = 16, cex = 0.95)
  points(rep(1650, 20), starting_scores2g, col = cols, pch = 16, cex = 0.95)
  points(rep(1750, 20), starting_scores2h, col = cols, pch = 16, cex = 0.95)
  points(rep(1850, 20), starting_scores2i, col = cols, pch = 16, cex = 0.95)
  points(rep(1950, 20), starting_scores2j, col = cols, pch = 16, cex = 0.95)
  points(rep(4000, 20), starting_scores3, col = cols, pch = 16, cex = 0.95)
  for (i in 1:ncol(X_sim)) {
    lines(1:1000, elo_sim[1:1000, i], col = cols[i], lty = lty_period1[i], lwd = 2)
    lines(1001:4000, elo_sim[1001:4000, i], col = cols[i], lty = lty_period2[i], lwd = 2)
  }
  points(rep(-139, 22), starting_scores1, col = cols_without_alpha, pch = rev(LETTERS[1:22])[order(starting_scores1)])
  points(rep(4140, 20), starting_scores3, col = cols_without_alpha, pch = rev(LETTERS[1:22])[order(starting_scores3)])
  ## Calculate recovery of true underlying 'hierarchy':
  elo_2nd_period <- est_hierarchy <- error <- error_corrected_for_sd <- elo_sim[1:4000, 1:20]
  ## 1:
  truth <- order(starting_scores1[1:20])
  for (i in 1:1000) {
    est_hierarchy[i, ] <- order(elo_2nd_period[i, ])
    error[i, ] <- error_corrected_for_sd[i, ] <- est_hierarchy[i, ] - truth
    error_corrected_for_sd[i, which(abs(error[i, ]) < 3)] <- 0
  }
  ## 2a:
  truth <- order(starting_scores2a)
  for (i in 1001:1100) {
    est_hierarchy[i, ] <- order(elo_2nd_period[i, ])
    error[i, ] <- error_corrected_for_sd[i, ] <- est_hierarchy[i, ] - truth
    error_corrected_for_sd[i, which(abs(error[i, ]) < 3)] <- 0
  }
  ## 2b:
  truth <- order(starting_scores2b)
  for (i in 1101:1200) {
    est_hierarchy[i, ] <- order(elo_2nd_period[i, ])
    error[i, ] <- error_corrected_for_sd[i, ] <- est_hierarchy[i, ] - truth
    error_corrected_for_sd[i, which(abs(error[i, ]) < 3)] <- 0
  }
  ## 2c:
  truth <- order(starting_scores2c)
  for (i in 1201:1300) {
    est_hierarchy[i, ] <- order(elo_2nd_period[i, ])
    error[i, ] <- error_corrected_for_sd[i, ] <- est_hierarchy[i, ] - truth
    error_corrected_for_sd[i, which(abs(error[i, ]) < 3)] <- 0
  }
  ## 2d:
  truth <- order(starting_scores2d)
  for (i in 1301:1400) {
    est_hierarchy[i, ] <- order(elo_2nd_period[i, ])
    error[i, ] <- error_corrected_for_sd[i, ] <- est_hierarchy[i, ] - truth
    error_corrected_for_sd[i, which(abs(error[i, ]) < 3)] <- 0
  }
  ## 2e:
  truth <- order(starting_scores2e)
  for (i in 1401:1500) {
    est_hierarchy[i, ] <- order(elo_2nd_period[i, ])
    error[i, ] <- error_corrected_for_sd[i, ] <- est_hierarchy[i, ] - truth
    error_corrected_for_sd[i, which(abs(error[i, ]) < 3)] <- 0
  }
  ## 2f:
  truth <- order(starting_scores2f)
  for (i in 1501:1600) {
    est_hierarchy[i, ] <- order(elo_2nd_period[i, ])
    error[i, ] <- error_corrected_for_sd[i, ] <- est_hierarchy[i, ] - truth
    error_corrected_for_sd[i, which(abs(error[i, ]) < 3)] <- 0
  }
  ## 2g:
  truth <- order(starting_scores2g)
  for (i in 1601:1700) {
    est_hierarchy[i, ] <- order(elo_2nd_period[i, ])
    error[i, ] <- error_corrected_for_sd[i, ] <- est_hierarchy[i, ] - truth
    error_corrected_for_sd[i, which(abs(error[i, ]) < 3)] <- 0
  }
  ## 2h:
  truth <- order(starting_scores2h)
  for (i in 1701:1800) {
    est_hierarchy[i, ] <- order(elo_2nd_period[i, ])
    error[i, ] <- error_corrected_for_sd[i, ] <- est_hierarchy[i, ] - truth
    error_corrected_for_sd[i, which(abs(error[i, ]) < 3)] <- 0
  }
  ## 2i:
  truth <- order(starting_scores2i)
  for (i in 1801:1900) {
    est_hierarchy[i, ] <- order(elo_2nd_period[i, ])
    error[i, ] <- error_corrected_for_sd[i, ] <- est_hierarchy[i, ] - truth
    error_corrected_for_sd[i, which(abs(error[i, ]) < 3)] <- 0
  }
  ## 2j:
  truth <- order(starting_scores2j)
  for (i in 1901:2000) {
    est_hierarchy[i, ] <- order(elo_2nd_period[i, ])
    error[i, ] <- error_corrected_for_sd[i, ] <- est_hierarchy[i, ] - truth
    error_corrected_for_sd[i, which(abs(error[i, ]) < 3)] <- 0
  }
  ## 3:
  truth <- order(starting_scores3)
  for (i in 2001:4000) {
    est_hierarchy[i, ] <- order(elo_2nd_period[i, ])
    error[i, ] <- error_corrected_for_sd[i, ] <- est_hierarchy[i, ] - truth
    error_corrected_for_sd[i, which(abs(error[i, ]) < 3)] <- 0
  }
  sum_abs_error <- apply(error, MAR = 1, FUN = function(x){sum(abs(x))})
  sum_abs_error_corrected_for_sd <- apply(error_corrected_for_sd, MAR = 1, FUN = function(x){sum(abs(x))})
  plot(1:4000, sum_abs_error/length(truth), type = "s", bty = "n", xlab = "Interaction index", 
       ylab = "Ranks MAE", yaxt = "n", lwd = 1, ylim = c(0, 10))
  axis(2, las = 2)
  # lines(1:4000, sum_abs_error_corrected_for_sd/length(truth), lty = 2, type = "s", col = rgb(0, 0.5, 0.5), lwd = 1)
}
dev.off()
shell.exec(paste0("simulation4_", today, ".pdf"))