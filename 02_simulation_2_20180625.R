## ###############################
## External Takeover Simulation ##
## Uncertainty plots            ##
## ###############################
## Simulate 2000 interactions (balanced) and show the uncertainty of the 
## starting scores by uncertainty intervals.

starting_scores1 <- seq(-6, 6, length = 10)
starting_scores2 <- c(starting_scores1, max(starting_scores1))
starting_scores2[which.max(starting_scores1)] <- min(starting_scores1) - diff(starting_scores1)[1]

cols <- viridis(n = length(starting_scores2), option = "A", begin = 0, end = 0.8, alpha = 0.5)
cols_without_alpha <- viridis(n = length(starting_scores2), option = "A", begin = 0, end = 0.8)

## Estimation:
chains <- 1; thin <- 1; iter <- 1500; warmup <- 500
iter_used <- iter - warmup
repeats <- 2
m_list <- X_sim_list <- vector("list", repeats)
period_length <- 2000
set.seed(1604)
for (rep_index in 1:repeats) {
  X_sim1 <- simulate_external_takeover_data(starting_scores = starting_scores1, n_interaction = period_length)
  X_sim2 <- simulate_external_takeover_data(starting_scores = starting_scores2, n_interaction = period_length)
  X_sim <- rbind(cbind(X_sim1, 0), X_sim2)
  X_sim_list[[rep_index]] <- X_sim
  Ai <- apply(X_sim, MAR = 1, FUN = function(x){which(x == 1)})
  Bi <- apply(X_sim, MAR = 1, FUN = function(x){which(x == -1)})
  fit_dat <- list(N = nrow(X_sim), K = ncol(X_sim), Ai = Ai, Bi = Bi, 
                  y = rep(1, nrow(X_sim)), diff_f = 1)
  fit_dat$presence <- 1 + 0 * X_sim
  fit_dat$presence[1:period_length, ncol(X_sim)] <- 0
  m_list[[rep_index]] <- stan(file = 'elo_score_USE_THIS.stan', data = fit_dat, 
                              iter = iter*thin, chains = chains, thin = thin, warmup = warmup*thin, 
                              control = list(adapt_delta = 0.95))
}
## Calculate credible intervals:
## Comment out (until including 'save(...)'-commands to save time (after first time having saved the results)):
Elo_pA <- function(EloStart_logk, X, show_k = FALSE, presence){
  EloStart <- EloStart_logk[-length(EloStart_logk)]
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
lower1_list <- lower2_list <- upper1_list <- upper2_list <-
  med_list <- elo_list_list <- vector("list", repeats)
for (rep_index in 1:repeats) {
  X_sim <- X_sim_list[[rep_index]]
  draws_starting_scores <- extract(m_list[[rep_index]])[[4]]
  draws_k <- extract(m_list[[rep_index]])[[2]]
  lower1_list[[rep_index]] <- upper1_list[[rep_index]] <-
    lower2_list[[rep_index]] <- upper2_list[[rep_index]] <-
    med_list[[rep_index]] <- vector("list", ncol(X_sim))
  for (j in 1:ncol(X_sim)) {
    cat(paste0(j, ":\n"))
    elo_list <- vector("list", nrow(draws_starting_scores))
    for (i in 1:nrow(draws_starting_scores)) {
      elo_list[[i]] <- Elo_pA(EloStart_logk = c(draws_starting_scores[i, ], log(draws_k[i])),
                              X = X_sim, presence = fit_dat$presence)$Elo
      ## progress trakler:
      cat(".")
      if ((i %% 80) == 0) {
        cat(".\n")
      }
    }
    cat("\n")
    aux <- NULL
    for (i in 1:length(elo_list)) {
      aux <- cbind(aux, elo_list[[i]][, j])
    }
    aux <- as.matrix(apply(aux, MAR = 1, FUN = quantile, probs = c(0.025, 0.1, 0.5, 0.9, 0.975)))
    med_list[[rep_index]][[j]] <- aux[3, ]
    upper1_list[[rep_index]][[j]] <- aux[1, ]
    upper2_list[[rep_index]][[j]] <- aux[2, ]
    lower2_list[[rep_index]][[j]] <- aux[4, ]
    lower1_list[[rep_index]][[j]] <- aux[5, ]
  }
  elo_list_list[[rep_index]] <- elo_list
  cat("\n")
}
save(med_list, file = paste0("sim2_med_list_", today, ".RData"))
save(lower1_list, file = paste0("sim2_lower1_list_", today, ".RData"))
save(lower2_list, file = paste0("sim2_lower2_list_", today, ".RData"))
save(upper1_list, file = paste0("sim2_upper1_list_", today, ".RData"))
save(upper2_list, file = paste0("sim2_upper2_list_", today, ".RData"))
save(elo_list_list, file = paste0("sim2_elo_list_list_", today, ".RData"))
## This will save time:
# load(file = paste0("sim2_med_list_", today, ".RData"))
# load(file = paste0("sim2_lower1_list_", today, ".RData"))
# load(file = paste0("sim2_lower2_list_", today, ".RData"))
# load(file = paste0("sim2_upper1_list_", today, ".RData"))
# load(file = paste0("sim2_upper2_list_", today, ".RData"))
# load(file = paste0("sim2_elo_list_list_", today, ".RData"))

## Plot results:
pdf(paste0("simulation2_", today, ".pdf"), height = 5, width = 7)
layout(mat = matrix(nrow = 2, ncol = 2, c(1, 3, 2, 4), byrow = T), 
       heights = c(0.65, 0.35), width = c(0.5, 0.5))
par(mar = c(4, 4, 1, 1))
for (rep_index in 1:repeats) {
  elo_sim <- do.call(what = cbind, med_list[[rep_index]])
  X_sim <- X_sim_list[[rep_index]]
  range_elo <- range(c(unlist(lower1_list[[rep_index]])), 
                     c(unlist(upper1_list[[rep_index]])))
  plot(elo_sim[, 1], type = "n", bty = "n", xlab = "", ylab = "", yaxt = "n", 
       ylim = range_elo, xaxt = "n", xlim = c(-100, 4100))
  draws_k <- extract(m_list[[rep_index]])[[2]]
  main_title <- paste("k: ", round(mean(draws_k), 2), ", with 95% CI: [", 
                      round(quantile(draws_k, probs = 0.025), 2), ", ", 
                      round(quantile(draws_k, probs = 0.975), 2), "]", sep = "")
  title(main = main_title)
  mtext(side = 1, text = "Interaction index", line = 2, cex = 1)
  mtext(side = 2, text = "Elo-rating", line = 2, cex = 1)
  axis(1, lab = FALSE, at = seq(0, 4000, by = 1000))
  axis(2, lab = FALSE, at = seq(-5, 5, by = 5))
  mtext(side = 1, at = seq(0, 4000, by = 1000), text = seq(0, 4000, by = 1000), 
        cex = 1, line = 0.7)
  mtext(side = 2, at = seq(-5, 5, by = 5), text = seq(-5, 5, by = 5), cex = 1, 
        line = 0.7, las = 2)
  lines(c(2000.5, 2000.5), range(elo_sim), lty = 2)
  points(rep(1, 10), starting_scores1, col = cols, pch = 16)
  points(rep(-139, 10), starting_scores1, col = cols_without_alpha, pch = rev(LETTERS[1:10])[order(starting_scores1)])
  points(rep(4000, 11), starting_scores2, col = cols, pch = 16)
  points(rep(4140, 11), starting_scores2, col = cols_without_alpha, pch = c(rev(LETTERS[2:10]), "A", "K"))
  pres_here <- 1:4000
  for (j in c(2, 6, 10)) {
    polygon(c(pres_here, rev(pres_here)), 
            c(lower1_list[[rep_index]][[j]], rev(upper1_list[[rep_index]][[j]])), 
            col = cols[j], border = NA)
    lines(pres_here, elo_sim[, j], col = cols_without_alpha[j], 
          lwd = 1, lty = 1)
  }
  j <- 11
  polygon(c(2001:4000, rev(2001:4000)), 
          c(lower1_list[[rep_index]][[j]][2001:4000], rev(upper1_list[[rep_index]][[j]][2001:4000])), 
          col = cols[j], border = NA)
  lines(2001:4000, elo_sim[2001:4000, j], col = cols_without_alpha[j], 
        lwd = 1, lty = 1)
  ## Calculate recovery of true underlying 'hierarchy':
  elo_both_periods <- est_hierarchy <- error <- elo_sim[1:4000, ]
  error[1:2000, 11] <- NA ## overwrite non-valid values with NA
  ## Period 1:
  truth <- order(starting_scores1) ## length: 10
  for (i in 1:2000) {
    est_hierarchy[i, 1:10] <- order(elo_both_periods[i, 1:10])
    error[i, 1:10] <- est_hierarchy[i, 1:10] - truth
  }
  ## Period 2:
  truth <- order(starting_scores2) ## length: 11
  for (i in 2001:4000) {
    est_hierarchy[i, ] <- order(elo_both_periods[i, ])
    error[i, ] <- est_hierarchy[i, ] - truth
  }
  sum_abs_error <- apply(error, MAR = 1, FUN = function(x){sum(abs(x), na.rm = T)})
  sum_abs_error[1:2000] <- sum_abs_error[1:2000]/10
  sum_abs_error[2001:4000] <- sum_abs_error[2001:4000]/11
  plot(1:4000, sum_abs_error, type = "s", bty = "n", ylim = c(0, 2), 
       xlab = "Interaction index", ylab = "Ranks MAE", yaxt = "n", lwd = 1)
  axis(2, las = 2)
}
dev.off()
shell.exec(paste0("simulation2_", today, ".pdf"))