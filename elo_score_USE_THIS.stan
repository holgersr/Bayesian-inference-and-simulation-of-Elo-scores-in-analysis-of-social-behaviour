functions {
  real[] ProbFunction(real[] EloStart, real k, matrix presence, int N, int K, int[] Ai, int[] Bi, real diff_f) {
    real result[N];
    real toAdd;
    //real aux_mean = 0.0;
    vector[K] EloNow;
    for (j in 1:K) {
      EloNow[j] = EloStart[j];
    }
    for (i in 1:N) {
      // centering:
      EloNow = EloNow - dot_product(row(presence,i),EloNow)/sum(row(presence,i));
      // likelihood contribution:
      result[i] = 1/(1 + exp(diff_f * (EloNow[Bi[i]] - EloNow[Ai[i]])));
      // update addend:
      toAdd = (1 - result[i]) * k;
      // update:
      EloNow[Ai[i]] = EloNow[Ai[i]] + toAdd;
      EloNow[Bi[i]] = EloNow[Bi[i]] - toAdd;
      }
    return result;
  }
} 
data {
  int<lower=1> N; // number of encounters
  int<lower=1> K; // number of individuals
  int<lower=1> Ai[N]; // winner's index
  int<lower=1> Bi[N]; // losers's index
  matrix[N, K] presence;
  int<lower=0> y[N]; // always 1
  real<lower=0> diff_f; // Elo Score difference factor
}
parameters {
  real EloStart_raw[K];
  real<lower=0.0> k_raw;
  real<lower=0.0> sigma_raw;
}
transformed parameters {
  real EloStart[K];
  real<lower=0.0> k;
  for (i in 1:K) {
      EloStart[i] = EloStart_raw[i] - mean(EloStart_raw);
    }
  for (i in 1:K) {
      EloStart[i] = EloStart[i]/diff_f;
    }
  k = k_raw/diff_f;
  }
model {
  k_raw ~ normal(0, 1);
  sigma_raw ~ normal(0, 1);
  EloStart_raw ~ normal(0, sigma_raw); 
  y ~ bernoulli(ProbFunction(EloStart, k, presence, N, K, Ai, Bi, diff_f));
  }
generated quantities{
  real<lower=0.0> sigma;
  sigma = sigma_raw/diff_f;
}
