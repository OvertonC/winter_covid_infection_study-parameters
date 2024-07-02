// distribution: lognormal

// truncation: right

// doubly interval censored

data {
  int<lower = 0> N;                           // number unique (A minus, A minus, B minus, B plus) combos
  vector<lower = 0>[N] a_minus;               // lower limit of event A
  vector<lower = 0>[N] a_plus;                // upper limit of event A
  vector<lower = 0>[N] b_minus;               // lower limit of event B
  vector<lower = 0>[N] b_plus;                // upper limit of event B
  vector<lower = 0>[3] incubation;            // [1] = T/F inclusion of time from O to A
                                              // (incubation period); [2] = mean [3] = sd
  real<lower = 0> upper_bound;                // the latest time of observation
  int A;                                      // number of clusters
  array[N] int<lower=1,upper=A> cluster;      // cluster ID variable
}

parameters {
  real log_mean_intercept;
  vector[A] beta_log_mean_cluster;
  real<lower=0> sigma_log_mean_cluster;
  real log_sd_intercept;
  vector[A] beta_log_sd_cluster;
  real<lower=0> sigma_log_sd_cluster;
  vector<lower = 0, upper = 1>[N] a_window;   // where time a lies in the event A window
  vector<lower = 0, upper = 1>[N] b_window;   // where time b lies in the event B window
  vector<lower = 0>[N] t0;                    // time from O to A
}

transformed parameters {
  vector[A] log_mean = log_mean_intercept + beta_log_mean_cluster * sigma_log_mean_cluster;
  vector[A] log_sd = log_sd_intercept + beta_log_sd_cluster * sigma_log_sd_cluster;
  vector[A] mu = 2*log(exp(log_mean)) - 0.5*log(exp(log_sd).^2 + exp(log_mean).^2);
  vector[A] sigma = sqrt(log(1 + (exp(log_sd).^2)./(exp(log_mean).^2) ));

  vector<lower = min(a_minus), upper = max(a_plus)>[N] a;
  vector<lower = min(b_minus), upper = max(b_plus)>[N] b;
  vector[N] ub;
  b = b_minus + (b_plus - b_minus) .* b_window;
  for (n in 1:N) {
    ub[n] = min([a_plus[n], b[n]]');
  }
  a = a_minus + (ub - a_minus) .* a_window;
}

model {
  log_mean_intercept ~ normal(1,2);
  log_sd_intercept ~ normal(1.6, 1);
  beta_log_mean_cluster ~ normal(0,1);
  sigma_log_mean_cluster ~ exponential(1);
  beta_log_sd_cluster ~ normal(0,1);
  sigma_log_sd_cluster ~ exponential(1);

  if (incubation[1] == 1) {
    t0 ~ lognormal(incubation[2], incubation[3]);
  } else {
    t0 ~ normal(0, 1e-10);
  }
  target += lognormal_lpdf((b - a + t0) | mu[cluster], sigma[cluster])
  - log(lognormal_cdf((upper_bound - a + t0) | mu[cluster], sigma[cluster]) - lognormal_cdf(1 | mu[cluster], sigma[cluster]));
  }

generated quantities {
  real<lower=0> mean_all_clusters_ = exp(log_mean_intercept);
  real<lower=0> sd_all_clusters_ = exp(log_sd_intercept);
  vector<lower = 0>[A] mean_ = exp(log_mean);
  vector<lower = 0>[A] sd_ = exp(log_sd);
  vector<lower = 0>[A] limit_val_ = exp(mu + sigma*(1.959963984540));
  vector[N] log_likelihood;
  for (n in 1:N) {
    log_likelihood[n] = lognormal_lpdf(b[n] - a[n] + t0[n] | mu[cluster[n]], sigma[cluster[n]])
    - log(lognormal_cdf((upper_bound - a[n] + t0[n]) | mu[cluster[n]], sigma[cluster[n]]) - lognormal_cdf(1 | mu[cluster[n]], sigma[cluster[n]]));
  }
}
