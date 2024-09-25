#include functions.stan

data {
  int<lower=1> N_total;
  int<lower=1> N_obs;
  int<lower=1> num_timeseries;

  int<lower=1> p; 
  int<lower=1> q; 
  array[N_total] int observed_TT;
  array[num_timeseries] int N_total_list;

  array[p+1, N_obs] int<lower=-1> Y_obs;

  array[N_total] vector[q] FF;
  array[N_total] matrix[q,q] GG;
  array[N_total] real gamma;
  
  array[num_timeseries] matrix[q,p] M0;
  array[num_timeseries] cov_matrix[q] C0;
  cov_matrix[p] Xi0;
  real<lower=p-1> upsilon0;

  real<lower=0> a;
  real<lower=0> b;
  real<lower=0> c;
  real<lower=0> d;
}

parameters {
  matrix[p, N_obs] eta;
  real<lower=0> W_theta;
  real<lower=0> W_alpha;  
}

transformed parameters {
  array[N_obs] simplex[p+1] pi;
  for (j in 1:N_obs) {
      pi[j] = softmax_id(eta[,j]);
  }
  array[N_total] matrix[q, q] WW;  
  for (n in 1:N_total) {
    WW[n] = identity_matrix(q);  
    WW[n][1,1] = W_theta;
    WW[n][2,2] = W_alpha;
  }
}

model {
  int pos = 1;

  W_theta ~ inv_gamma(a,b);
  W_alpha ~ inv_gamma(c,d);
  target += gmdlm_lpdf(eta | FF, GG, WW, gamma, M0, C0, Xi0, upsilon0, observed_TT, N_total_list, num_timeseries);
  for (i in 1:N_total) {
    if (observed_TT[i] == 1) {
      Y_obs[, pos] ~ multinomial(pi[pos]);
      pos = pos + 1;
    }
  }
}

generated quantities {
  array[N_total] matrix[q, p] theta;
  theta = gmdlm_smoothing_rng(eta, FF, GG, WW, gamma, M0, C0, Xi0, upsilon0, observed_TT, N_total_list, num_timeseries);
}
