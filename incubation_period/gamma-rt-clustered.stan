// distribution: gamma
// truncation:   right
// doubly interval censored
data {
 int<lower = 0> N;                     // number of records
 vector<lower = 0>[N] a_minus;         // lower limit of event A
 vector<lower = 0>[N] a_plus;          // upper limit of event A
 vector<lower = 0>[N] b_minus;         // lower limit of event B
 vector<lower = 0>[N] b_plus;          // upper limit of event B
 vector<lower = 0>[3] incubation;      // [1] = T/F inclusion of time from O to A
                                       //  (incubation period); [2] = mean [3] = sd
 real<lower = 0> upper_bound;          // the latest time of observation
 int A;                                // number of clusters
 array[N] int<lower=1,upper=A> cluster;// cluster ID variable
}

parameters {
 real log_mean_intercept;                                 // the distribution mean
 vector[A] beta_log_mean_cluster; //per-cluster random effect on log mean
 real<lower=0> sigma_log_mean_cluster;
 real log_sd_intercept;  // log of distribution standard deviation
 vector[A] beta_log_sd_cluster; //per-cluster random effect on log sd
 real<lower = 0> sigma_log_sd_cluster;
 vector<lower = 0, upper = 1>[N] a_window;  // where time a lies in the event A window
 vector<lower = 0, upper = 1>[N] b_window;  // where time b lies in the event B window
 vector<lower = 0>[N] t0;                             // time from O to A
}

transformed parameters {
 vector[A] logmean_ = log_mean_intercept + beta_log_mean_cluster*sigma_log_mean_cluster;
 vector<lower = 0>[A] sd_ = exp(log_sd_intercept + beta_log_sd_cluster*sigma_log_sd_cluster);
 vector<lower = 0>[A] mean_ = exp(logmean_);
 vector<lower = 0>[A] alpha = (mean_./sd_)^2;
 vector<lower = 0>[A] beta = mean_./(sd_.^2);
 vector<lower = min(a_minus), upper = max(a_plus)>[N] a;
 vector<lower = min(b_minus), upper = max(b_plus)>[N] b;
 vector[N] ub;
 b = b_minus + (b_plus - b_minus) .* b_window;
 for (n in 1:N)
   ub[n] = min([a_plus[n], b[n]]');
 a = a_minus + (ub - a_minus) .* a_window;

}

model {
 log_mean_intercept ~ normal(1,2);
 log_sd_intercept ~ normal(1.6, 1);
 beta_log_mean_cluster ~ std_normal();
 sigma_log_mean_cluster ~ exponential(1);
 beta_log_sd_cluster ~ std_normal();
 sigma_log_sd_cluster ~ exponential(1);

 if (incubation[1] == 1) {
   t0 ~ lognormal(incubation[2], incubation[3]);
 } else {
   t0 ~ normal(0, 1e-10);
 }
 target += gamma_lpdf((b - a + t0) | alpha[cluster], beta[cluster])
 - log(gamma_cdf((upper_bound - a + t0) | alpha[cluster], beta[cluster]) - gamma_cdf(1 | alpha[cluster], beta[cluster]));
}

generated quantities {
 vector<lower=0>[A] sd2_;
 sd2_ = sqrt(alpha./beta./beta);
 vector[N] log_likelihood;
 for (n in 1:N) {
   log_likelihood[n] = gamma_lpdf(b[n] - a[n] + t0[n] | alpha[cluster[n]], beta[cluster[n]])
   - log(gamma_cdf((upper_bound - a[n] + t0[n]) | alpha[cluster[n]], beta[cluster[n]]) - gamma_cdf(1 | alpha[cluster[n]], beta[cluster[n]]));
  }
}
