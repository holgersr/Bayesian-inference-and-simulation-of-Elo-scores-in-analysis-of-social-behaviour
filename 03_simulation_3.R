## ####################
## Coalitionary Leap ##
## ####################
## 4 groups of 5 individuals
## March 19, 2017

cols <- c(viridis(n = 5, option = "A", alpha = 0.5, begin = 0.15-0.1, end = 0.2-0.1), 
          viridis(n = 5, option = "A", alpha = 0.5, begin = 0.35-0.1, end = 0.4-0.1), 
          viridis(n = 5, option = "A", alpha = 0.5, begin = 0.55-0.1, end = 0.6-0.1), 
          viridis(n = 5, option = "A", alpha = 0.5, begin = 0.75-0.1, end = 0.8-0.1), 
          viridis(n = 5, option = "A", alpha = 0.5, begin = 0.95-0.1, end = 1-0.1))

x <- 2*seq(6, -6, length = 25)
starting_scores1 <- x
starting_scores2 <- x[c(6:15, 1:5, 16:25)]
chains <- 1; thin <- 10; iter <- 1500; warmup <- 500
iter_used <- iter - warmup

setwd(path_here)

repeats <- 2

set.seed(1604)
x <- seq(6, -6, length = 25)
fit_list <- X_sim_list <- vector("list", repeats)
for(rep_index in 1:repeats){
  starting_scores1 <- x 
  starting_scores2 <- x[c(6:15, 1:5, 16:25)]
  X_sim1 <- simulate_coalitionary_leap_data(starting_scores = starting_scores1, 
                                            n_interaction = 1000)
  X_sim2 <- simulate_coalitionary_leap_data(starting_scores = starting_scores2, 
                                            n_interaction = 2000)
  X_sim <- rbind(X_sim1, X_sim2)
  X_sim_list[[rep_index]] <- X_sim
  
  Ai <- apply(X_sim, MAR = 1, FUN = function(x){which(x == 1)})
  Bi <- apply(X_sim, MAR = 1, FUN = function(x){which(x == -1)})
  fit_dat <- list(N = nrow(X_sim), K = ncol(X_sim), Ai = Ai, Bi = Bi, 
                  y = rep(1, nrow(X_sim)), diff_f = 1)
  fit_dat$presence <- 1 + 0 * X_sim
  fit_list[[rep_index]] <- stan(file = 'elo_score_USE_THIS.stan', data = fit_dat, 
                                iter = iter*thin, chains = chains, 
                                thin = thin, warmup = warmup*thin, 
                                control = list(adapt_delta = 0.95))
}
## ################
## Plot results: ##
## ################
Elo_pA <- function(EloStart_logk, X, show_k = FALSE, presence){
  EloStart <- EloStart_logk[-length(EloStart_logk)]
  k <- exp(EloStart_logk[length(EloStart_logk)])
  if(show_k){
    cat(paste0(round(k, 3), paste(rep(" ", 20), collapse = " ")))
    cat("\r")
  }
  EloNow <- EloStart
  Elo <- matrix(nrow = nrow(X), ncol = length(EloStart), 0)
  colnames(Elo) <- colnames(X)
  pA <- rep(0, nrow(X))
  for(i in 1:nrow(X)){
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

pdf("simulation3.pdf", height = 5, width = 7)
layout(mat = matrix(nrow = 2, ncol = 2, c(1, 3, 2, 4), byrow = T), 
       heights = c(0.65, 0.35), width = c(0.5, 0.5))
par(mar = c(4, 4, 1, 1))
for(rep_index in 1:repeats){
  bayes_model3 <- fit_list[[rep_index]]
  X_sim <- X_sim_list[[rep_index]]
  bayes_starting_scores <- summary(bayes_model3)$summary[1:ncol(X_sim)]
  bayes_k <- summary(bayes_model3)$summary[ncol(X_sim)+1]
  elo_sim <- Elo_pA(EloStart_logk = c(bayes_starting_scores, log(bayes_k)), 
                    X = X_sim, presence = 1 + 0 * X_sim)$Elo
  plot(elo_sim[, 1], type = "n", ylim = range(c(elo_sim, starting_scores1, 
                                                starting_scores2)), 
       bty = "n", yaxt = "n", xlab = "Interaction index", ylab = "Elo score")
  main_title <- paste("k: ", round(bayes_k, 2), ", with 95% CI: [", 
                      paste0(unname(round(summary(bayes_model3)$summary[ncol(X_sim)+1, c(4, 8)], 2)), 
                             collapse = ", "), "]", sep = "")
  title(main = main_title)
  axis(2, lab = FALSE, at = seq(-5, 5, by = 5))
  mtext(side = 2, at = seq(-5, 5, by = 5), text = seq(-5, 5, by = 5), cex = 1, 
        line = 0.7, las = 2)
  lty_period1 <- c(rep(1, 25))
  lty_period2 <- c(rep(1, 25))
  lwd_period1 <- lwd_period2 <- c(rep(2, 25))
  points(rep(1, 25), starting_scores1, col = cols, pch = 16)
  points(rep(3000, 25), starting_scores2, col = cols, pch = 16)
  for(i in 1:ncol(X_sim)){
    lines(1:1001, elo_sim[1:1001, i], col = cols[i], lty = lty_period1[i], 
          lwd = lwd_period1[i])
    lines(1001:3000, elo_sim[1001:3000, i], col = cols[i], lty = lty_period2[i], 
          lwd = lwd_period2[i])
  }
  lines(rep(1001, 2), range(starting_scores1), lty = 2)
  ## Calculate recovery of true underlying 'hierarchy':
  elo_2nd_period <- est_hierarchy <- error <- 
    error_corrected_for_sd <- elo_sim[1:3000, ]
  truth <- order(starting_scores1)
  for(i in 1:1000){
    est_hierarchy[i, ] <- order(elo_2nd_period[i, ])
    error[i, ] <- error_corrected_for_sd[i, ] <- est_hierarchy[i, ] - truth
    error_corrected_for_sd[i, which(abs(error[i, ]) < 3)] <- 0
  }
  truth <- order(starting_scores2)
  for(i in 1001:3000){
    est_hierarchy[i, ] <- order(elo_2nd_period[i, ])
    error[i, ] <- error_corrected_for_sd[i, ] <- est_hierarchy[i, ] - truth
    error_corrected_for_sd[i, which(abs(error[i, ]) < 3)] <- 0
  }
  sum_abs_error <- apply(error, MAR = 1, FUN = function(x){sum(abs(x))})
  sum_abs_error <- sum_abs_error/length(truth)
  plot(1:3000, sum_abs_error, type = "s", bty = "n", xlab = "Interaction index", 
       ylab = "Ranks MAE", yaxt = "n", lwd = 1, ylim = c(0, max(sum_abs_error)))
  axis(2, las = 2)
  }
dev.off()
