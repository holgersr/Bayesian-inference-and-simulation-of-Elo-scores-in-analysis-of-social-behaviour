## ############################################################################
## Re-analysis of eastern chimpanzee (Pan troglodytes schweinfurthii)        ##
## agonistic interaction data (female subgroup),                             ##
## published under: http://datadryad.org/resource/doi:10.5061/dryad.r4g74    ##
## Holger Sennhenn-Reulen, Adeelia Goffe, Julia Fischer                      ##
## June 22, 2017                                                             ##
## ############################################################################
library("rstan")
library("viridis")
library("loo")
library("ggplot2")
library("plyr")

source("00a_functions_ml_approach.R")

## Female agonistic interaction data:
d <- read.csv(file = "foerster2016_female_ago.csv", header = T, sep = ";", 
              stringsAsFactors = FALSE)
pres <- read.csv("foerster2016_female_presence.csv", header = T, sep = ";", 
                 stringsAsFactors = FALSE)
length(unique(c(d$Winner, d$Loser))) ## 44 individuals were observed.
X <- matrix(nrow = nrow(d), ncol = length(unique(c(d$Winner, d$Loser))), 0)
for(i in 1:nrow(d)){
  X[i, d$Winner[i]] <- 1
  X[i, d$Loser[i]] <- -1
}
apply(X, MAR = 2, FUN = function(x){sum(x == 1)})[9]
apply(X, MAR = 2, FUN = function(x){sum(x == -1)})[9]

X <- X[101:nrow(X), ] ## First 100 interaction were removed by Foerster et al., 2016

## #######################################################
## Foerster 2016 approch (direct log. Lik. maximation): ##
## #######################################################
## For convenience, results were only calculated once, and loaded thereafter.
## foerster_result_model3 <- optim(par = c(rep(0, ncol(X)), log(1)), fn = logLik, 
##                                 X = X, show_k = TRUE, method='BFGS', hessian = TRUE, 
##                                 control = list(maxit = 10000, reltol = 1e-10))
## foerster_result_model3_001factor <- optim(par = c(rep(0, ncol(X)), log(100)), 
##                                           fn = logLik_001factor, X = X, show_k = TRUE, 
##                                           method='BFGS', hessian = TRUE, 
##                                           control = list(maxit = 10000, reltol = 1e-10))
## save(foerster_result_model3, file = "foerster_result_model3.RData")
## save(foerster_result_model3_001factor, file = "foerster_result_model3_001factor.RData")
load(file = "foerster_result_model3.RData")
load(file = "foerster_result_model3_001factor.RData")

(foerster_result_model3_coef <- foerster_result_model3$par[-length(foerster_result_model3$par)])
(k_model3_foerster <- exp(foerster_result_model3$par[length(foerster_result_model3$par)]))
(foerster_result_model3_001factor_coef <- foerster_result_model3_001factor$par[-length(foerster_result_model3_001factor$par)])
(k_model3_001factor_foerster <- exp(foerster_result_model3_001factor$par[length(foerster_result_model3_001factor$par)]))

## ################################
## Presence matrix for centering ##
## ################################
## Manage observation and presence dates:
presence <- read.csv("../data/foerster2016_female_presence.csv", sep = ";")
presence_dates <- as.Date(presence$Date, format = "%d.%m.%Y")
presence <- presence[, -1]
observation_dates <- as.Date(d$Date[101:nrow(d)], format = "%d.%m.%Y")
aux_index <- which(presence_dates %in% observation_dates)
presence_dates <- presence_dates[aux_index]
presence <- presence[aux_index, ]

table(sort(unique(presence_dates)) == sort(unique(observation_dates)))
presence_dates_index_in_observation_dates <- rep(0, length(observation_dates))
for(i in 1:length(observation_dates)){
  aux_index <- which(presence_dates == observation_dates[i])
  presence_dates_index_in_observation_dates[i] <- aux_index
}
presence_for_estimation <- presence[presence_dates_index_in_observation_dates, ]
rm(aux_index)

## ###############################
## Bayesian inference approach: ##
## ###############################
Ai <- apply(X, MAR = 1, FUN = function(x){which(x == 1)})
Bi <- apply(X, MAR = 1, FUN = function(x){which(x == -1)})
chains <- 16; iter <- 500; warmup <- 250; thin <- 1
## Data provided to the STAN call:
fit_dat <- list(N = nrow(X), K = ncol(X), Ai = Ai, Bi = Bi, y = rep(1, nrow(X)), 
                diff_f = NULL)
## STAN call:
fit_dat$diff_f <- 1 ## Elo-score difference factor
fit_dat$presence <- as.matrix(presence_for_estimation)
fit <- stan(file = 'elo_score_USE_THIS.stan', data = fit_dat,
            iter = iter*thin, chains = chains, thin = thin, warmup = warmup*thin, 
            control = list(adapt_delta = 0.95), seed = 123)
summary(fit)$summary

pdf("prior_posterior_k_sigma.pdf", height = 5, width = 7) 
par(mfrow = c(1, 2))
draws <- extract(fit)[["k"]]
hist(draws, freq = F, yaxt = "n", col = "grey", border = "grey", main = "", 
     xlab = "k")
axis(2, las = 2)
x <- seq(min(draws), max(draws), length = 100)
lines(x, 2*dnorm(x))
draws <- extract(fit)[["sigma"]]
hist(draws, freq = F, yaxt = "n", col = "grey", border = "grey", main = "", 
     xlab = expression(paste(sigma)))
axis(2, las = 2)
x <- seq(min(draws), max(draws), length = 100)
lines(x, 2*dnorm(x, sd = 1))
dev.off()

## #######################
## Elo score path plot: ##
## #######################
EloStart <- extract(fit)[['EloStart']]
k <- extract(fit)[['k']]

## work correctly with equal observation dates:
D <- observation_dates
aux_rle <- rle(as.numeric(D))
aux_rle$lengths
aux_index <- which(aux_rle$lengths > 1.5)
daytime <- rep("12:00", length(D))
for(i in aux_index){
  where <- cumsum(aux_rle$lengths)[i-1] + 1:aux_rle$lengths[i]
  aux_daytime <- seq(0, 24*60, length = aux_rle$lengths[i]+2)
  aux_daytime <- aux_daytime[-1]
  aux_daytime <- aux_daytime[-length(aux_daytime)]
  h <- floor(aux_daytime/60)
  m <- floor(aux_daytime %% 60)
  h <- formatC(h, width = 2, format = "d", flag = "0")
  m <- formatC(m, width = 2, format = "d", flag = "0")
  daytime[where] <- paste(h, m, sep = ":")
}
D <- as.POSIXct(gsub(paste(D, daytime, sep = " "), pattern = "-", replacement = ""), 
                format = "%Y%m%d %H:%M")
rm(where, h, m, aux_daytime, aux_rle, aux_index,daytime)
## Two functions needed for post-estimation elo -score calculation:
post_estimation_elo_score_calculation_a <- function(EloStart, k, N, K, Ai, Bi, diff_f, 
                                                         presence_for_estimation){
  cat("Calculate Elo scores for all posterior samples:\n")
  result_list <- group_mean_list <- vector("list", length(k))
  n_present <- unname(apply(presence_for_estimation, MAR = 1, FUN = sum))
  present_here <- apply(presence_for_estimation, MAR = 1, FUN = function(x){unname(which(x == 1))})
  for(iteration in 1:length(k)){
    result <- matrix(nrow = N + 1, ncol = K, NA)
    group_mean <- rep(0, N)
    aux <- NULL
    for(j in 1:K){
      result[1, j] <- EloStart[iteration, j]
    }
    for(i in 2:(N+1)){
      ## update addend:
      aux <- 1/(1 + exp(diff_f * (result[i-1, Bi[i]] - result[i-1, Ai[i]])))
      aux <- (1 - aux) * k[iteration]
      ## centering:
      group_mean[i-1] <- sum(result[i-1, present_here[[i-1]]])/n_present[i-1]
      result[i-1, ] <- result[i-1, ] - group_mean[i-1]
      for(j in 1:K){
        result[i, j] <- result[i-1, j]
      }
      ## update:
      result[i, Ai[i]] <- result[i, Ai[i]] + aux
      result[i, Bi[i]] <- result[i, Bi[i]] - aux
    }
    result <- result[-1, ]
    result_list[[iteration]] <- result
    group_mean_list[[iteration]] <- group_mean
    ## Progress tracker:
    cat("\r"); cat("Done:", round(100*iteration/length(k), 2), "%.                "); cat("\r")
  }
  return(list(result_list = result_list, group_mean_list = group_mean_list))
}
post_estimation_elo_score_calculation_b <- function(Elo_list, presence_for_estimation){
  cat("Calculate Elo score mean and quantiles per ID and date:\n")
  K <- ncol(presence_for_estimation)
  N <- nrow(presence_for_estimation)
  when <- id <- NULL
  for(j in 1:K){
    when_was_j_present <- which(presence_for_estimation[, j] == 1)
    id <- c(id, rep(j, length(when_was_j_present)))
    when <- c(when, when_was_j_present)
  }
  result <- data.frame(id = id,
                       Date = when,
                       mean_elo = NA,
                       q025_elo = NA,
                       q1_elo = NA,
                       q9_elo = NA,
                       q975_elo = NA, 
                       stringsAsFactors = FALSE)
  for(i in 1:nrow(result)){
    j <- result$id[i]
    day <- result$Date[i]
    Elo <- rep(NA, length(Elo_list))
    for(iteration in 1:length(Elo_list)){
      Elo[iteration] <- Elo_list[[iteration]][day, j]
    }
    result$mean_elo[i] <- mean(Elo)
    result$q025_elo[i] <- quantile(Elo, probs = 0.025)
    result$q1_elo[i] <- quantile(Elo, probs = 0.1)
    result$q9_elo[i] <- quantile(Elo, probs = 0.9)
    result$q975_elo[i] <- quantile(Elo, probs = 0.975)
    cat("\r"); cat("Done:", round(100*i/nrow(result), 1), "% (ID:", j, ")                     "); cat("\r")
  }
  return(result)
}

set.seed(123)
iteration_sample <- sample(1:length(k))[1:4000]
elo_a <- post_estimation_elo_score_calculation_a(EloStart = EloStart[iteration_sample, ], 
                                                      k = k[iteration_sample], 
                                                      N = fit_dat$N, K = fit_dat$K, 
                                                      Ai = fit_dat$Ai, Bi = fit_dat$Bi, 
                                                      diff_f = fit_dat$diff_f, 
                                                      presence_for_estimation = presence_for_estimation)

group_mean <- apply(do.call(cbind, elo_a[[2]]), MAR = 1, FUN = mean)
elo <- post_estimation_elo_score_calculation_b(Elo_list = elo_a[[1]], 
                                                    presence_for_estimation = presence_for_estimation)
rm(elo_a)
elo$id <- formatC(elo$id, width = max(nchar(as.character(elo$id))), 
                  format = "d", flag = "0")
Groupmean <- data.frame(Date = D, 
                        elo = cumsum(group_mean), 
                        stringsAsFactors = FALSE)
Groupmean <- ddply(Groupmean, c("Date"), summarize,  
                  mean_elo = mean(elo))
Groupmean$id_removed <- 1
Groupmean$id_removed <- as.factor(Groupmean$id_removed)

elo$Date <- D[elo$Date]
elo$Date <- as.Date(elo$Date)
Groupmean$Date <- as.Date(Groupmean$Date)

elo_copy <- elo
names(elo_copy)[1] <- "id_removed"
elo_poly <- NULL
for(i in unique(elo$id)){
  here <- subset(elo, id == i)
  elo_poly <- rbind(elo_poly, 
                    data.frame(Date = c(here$Date, rev(here$Date)),
                               outer_interval = c(here$q025_elo, rev(here$q975_elo)),
                               inner_interval = c(here$q1_elo, rev(here$q9_elo)),
                               id = i, stringsAsFactors = F))
}
foerster_results_matrix <- presence_for_estimation
for(j in 1:ncol(foerster_results_matrix)){
  foerster_results_matrix[, j] <- foerster_results_matrix[, j] * 
    foerster_result_model3_001factor_coef[j]/100
}
for(i in 1:nrow(foerster_results_matrix)){
  aux_index <- which(presence_for_estimation[i, ] == 1)
  aux_mean <- as.numeric(foerster_results_matrix[i, ])
  aux_mean <- aux_mean[aux_index]
  foerster_results_matrix[i, ] <- foerster_results_matrix[i, ] - mean(aux_mean)
}
for(i in 1:nrow(foerster_results_matrix)){
  for(j in 1:ncol(foerster_results_matrix)){
    if(presence_for_estimation[i, j] == 0){
      foerster_results_matrix[i, j] <- NA
    }
  }
}
foerster_et_al_results <- data.frame(elo = as.numeric(as.matrix(foerster_results_matrix)),
                                     id = rep(sort(unique(elo$id)), 
                                              each = nrow(foerster_results_matrix)),
                                     Date = as.Date(rep(D, ncol(foerster_results_matrix))))
foerster_et_al_results <- foerster_et_al_results[-which(is.na(foerster_et_al_results$elo)), ]
dd <- ddply(elo, c("id"), summarize, min_date = min(Date), max_date = max(Date))

col_foerster <- viridis(n = 1, option = "A", begin = 0.1, end = 0.1)
pb <- ggplot(elo, aes(x = Date, y = mean_elo, group = id)) +
  geom_line(data = Groupmean, aes(x = Date, y = mean_elo, group = id_removed), 
            size = 0.5, color = "black", linetype = "dotted") + 
  geom_line(data = foerster_et_al_results, aes(x = Date, y = elo, group = id), 
            size = 0.5, linetype = "solid", color = col_foerster) +
  geom_line(data = elo_copy, aes(x = Date, y = mean_elo, group = id_removed), 
            size = 0.2, linetype = "solid", color="grey") +  
  geom_polygon(data=elo_poly, mapping=aes(x = Date, y = outer_interval, 
                                          group = id, fill = factor(id))) + 
  geom_polygon(data=elo_poly, mapping=aes(x = Date, y = inner_interval, 
                                          group = id, fill = factor(id))) + 
  scale_fill_manual(values = rev(viridis(n = 44, option = "A", alpha = 0.5, 
                                         begin = 0.3, end = 0.8)), guide = FALSE) +
  geom_line(aes(color = factor(id)), size=0.5, alpha = 1) +
  scale_color_manual(values = rev(viridis(n = 44, option = "A", begin = 0.3, 
                                          end = 0.8)), guide = FALSE) + 
  ylab("") +
  xlab("") +
  theme_bw() +
  theme(panel.border=element_blank(), panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "grey", size = 0.1),
        strip.background=element_blank(), 
        axis.line = element_line(colour = "black")) + 
  facet_wrap(~ id) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  scale_y_continuous(breaks=seq(-5, 5, by = 5))
pb
## ggsave("foerster_female_data_elo_paths.pdf", pb, width=10, height=8)
## shell.exec("foerster_female_data_elo_paths.pdf")