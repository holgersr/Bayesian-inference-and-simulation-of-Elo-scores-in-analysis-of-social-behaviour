## Functions needed to implement the ML approach:
pAfun <- function(EloA, EloB){
  1/(1 + exp(EloB - EloA))
}
Elo_pA <- function(EloStart_logk, X, show_k = FALSE){
  EloStart <- EloStart_logk[-length(EloStart_logk)]
  EloStart <- EloStart - mean(EloStart)
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
    pA[i] <- pAfun(EloA = EloNow[A], EloB = EloNow[B])
    toAdd <- (1 - pA[i]) * k
    EloNow[A] <- EloNow[A] + toAdd
    EloNow[B] <- EloNow[B] - toAdd
    Elo[i, ] <- EloNow
  }
  return(list(pA = pA, Elo = Elo))
}
logLik <- function(EloStart_logk, X, show_k = FALSE){
  pA <- Elo_pA(EloStart_logk = EloStart_logk, X = X, show_k = show_k)$pA
  return(-sum(log(pA)))
}
logLik_model1 <- function(logk, X, show_k = FALSE){
  pA <- Elo_pA(EloStart_logk = c(rep(0, ncol(X)), logk), X = X, show_k = show_k)$pA
  return(-sum(log(pA)))
}
pAfun_001factor <- function(EloA, EloB){
  1/(1 + exp(0.01*(EloB - EloA)))
}
Elo_pA_001factor <- function(EloStart_logk, X, show_k = FALSE){
  EloStart <- EloStart_logk[-length(EloStart_logk)]
  EloStart <- EloStart - mean(EloStart)
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
    pA[i] <- pAfun_001factor(EloA = EloNow[A], EloB = EloNow[B])
    toAdd <- (1 - pA[i]) * k
    EloNow[A] <- EloNow[A] + toAdd
    EloNow[B] <- EloNow[B] - toAdd
    Elo[i, ] <- EloNow
  }
  return(list(pA = pA, Elo = Elo))
}
logLik_001factor <- function(EloStart_logk, X, show_k = FALSE){
  pA <- Elo_pA_001factor(EloStart_logk = EloStart_logk, X = X, show_k = show_k)$pA
  return(-sum(log(pA)))
}