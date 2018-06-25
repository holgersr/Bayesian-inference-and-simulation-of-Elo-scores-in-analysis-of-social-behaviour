## ##########################################################################
## Functions needed to simulate data according to the simulation scenarios ##
## described in the manuscript                                             ##
## ##########################################################################
simulate_foerster_data <- function(starting_scores, n_interaction = 200){
  n <- length(starting_scores)
  result <- data.frame(interaction_index = 1:n_interaction,
                       winner = rep(0, n_interaction),
                       loser = rep(0, n_interaction))
  for (i in 1:n_interaction) {
    id1 <- sample(1:n, size = 1)
    id2 <- sample(c(1:n)[-id1], size = 1, prob = 1/plogis(abs(starting_scores[-id1] - starting_scores[id1])))
    
    EloDiff = starting_scores[id2] - starting_scores[id1]
    p_id1 = 1/(1 + exp(EloDiff))
    
    id2_winner <- rbinom(size = 1, n = 1, prob = 1 - p_id1)
    
    result$winner[i] <- c(id1, id2)[id2_winner + 1]
    result$loser[i] <- c(id2, id1)[id2_winner + 1]
  }
  X <- matrix(nrow = n_interaction, ncol = n, 0)
  for (i in 1:n_interaction) {
    X[i, result$winner[i]] <- 1
    X[i, result$loser[i]] <- -1
  }
  return(X)
}

simulate_external_takeover_data <- function(starting_scores, n_interaction = 200){
  n <- length(starting_scores)
  result <- data.frame(interaction_index = 1:n_interaction,
                       winner = rep(0, n_interaction),
                       loser = rep(0, n_interaction))
  for (i in 1:n_interaction) {
    id1 <- sample(1:n, size = 1)
    id2 <- sample(c(1:n)[-id1], size = 1, prob = 1/plogis(abs(starting_scores[-id1] - starting_scores[id1])))
    EloDiff = starting_scores[id2] - starting_scores[id1]
    p_id1 = 1/(1 + exp(EloDiff))
    id2_winner <- rbinom(size = 1, n = 1, prob = 1 - p_id1)
    result$winner[i] <- c(id1, id2)[id2_winner + 1]
    result$loser[i] <- c(id2, id1)[id2_winner + 1]
  }
  X <- matrix(nrow = n_interaction, ncol = n, 0)
  for (i in 1:n_interaction) {
    X[i, result$winner[i]] <- 1
    X[i, result$loser[i]] <- -1
  }
  return(X)
}

simulate_coalitionary_leap_data <- function(starting_scores, n_interaction = 1000){
  n <- length(starting_scores)
  result <- data.frame(interaction_index = 1:n_interaction,
                       winner = rep(0, n_interaction),
                       loser = rep(0, n_interaction))
  for (i in 1:n_interaction) {
    id1 <- sample(1:n, size = 1)
    id2 <- sample(c(1:n)[-id1], size = 1, prob = 1/plogis(abs(starting_scores[-id1] - starting_scores[id1])))
    
    EloDiff = starting_scores[id2] - starting_scores[id1]
    p_id1 = 1/(1 + exp(EloDiff))
    
    id2_winner <- rbinom(size = 1, n = 1, prob = 1 - p_id1)
    
    result$winner[i] <- c(id1, id2)[id2_winner + 1]
    result$loser[i] <- c(id2, id1)[id2_winner + 1]
  }
  X <- matrix(nrow = n_interaction, ncol = n, 0)
  for (i in 1:n_interaction) {
    X[i, result$winner[i]] <- 1
    X[i, result$loser[i]] <- -1
  }
  return(X)
}

simulate_mortality_instability_recovery_data <- function(starting_scores, n_interaction = 200){
  n <- length(starting_scores)
  result <- data.frame(interaction_index = 1:n_interaction,
                       winner = rep(0, n_interaction),
                       loser = rep(0, n_interaction))
  for (i in 1:n_interaction) {
    id1 <- sample(1:n, size = 1)
    id2 <- sample(c(1:n)[-id1], size = 1, prob = 1/plogis(abs(starting_scores[-id1] - starting_scores[id1])))
    EloDiff = starting_scores[id2] - starting_scores[id1]
    p_id1 = 1/(1 + exp(EloDiff))
    id2_winner <- rbinom(size = 1, n = 1, prob = 1 - p_id1)
    result$winner[i] <- c(id1, id2)[id2_winner + 1]
    result$loser[i] <- c(id2, id1)[id2_winner + 1]
  }
  X <- matrix(nrow = n_interaction, ncol = n, 0)
  for (i in 1:n_interaction) {
    X[i, result$winner[i]] <- 1
    X[i, result$loser[i]] <- -1
  }
  return(X)
}

simulate_unbalanced_data <- function(strength, n_interaction = 200){
  n <- length(strength)
  result <- data.frame(interaction_index = 1:n_interaction,
                       winner = rep(0, n_interaction),
                       loser = rep(0, n_interaction))
  for (i in 1:n_interaction) {
    id1 <- sample(1:n, size = 1, prob = plogis(strength))
    id2 <- sample(c(1:n)[-id1], size = 1, prob = 1/plogis(abs(strength[-id1] - strength[id1])))
    EloDiff = strength[id2] - strength[id1]
    p_id1 = 1/(1 + exp(EloDiff))
    id2_winner <- rbinom(size = 1, n = 1, prob = 1 - p_id1)
    result$winner[i] <- c(id1, id2)[id2_winner + 1]
    result$loser[i] <- c(id2, id1)[id2_winner + 1]
  }
  X <- matrix(nrow = n_interaction, ncol = n, 0)
  for (i in 1:n_interaction) {
    X[i, result$winner[i]] <- 1
    X[i, result$loser[i]] <- -1
  }
  return(X)
}