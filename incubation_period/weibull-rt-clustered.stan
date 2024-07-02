// distribution: weibull
// truncation: right
// doubly interval censored
data {
 int<lower = 1> N;                           // number unique (A minus, A minus, B minus, B plus) combos
 vector<lower = 0>[N] a_minus;               // lower limit of event A
 vector<lower = 0>[N] a_plus;                // upper limit of event A
 vector<lower = 0>[N] b_minus;               // lower limit of event B
 vector<lower = 0>[N] b_plus;                // upper limit of event B
 vector<lower = 0>[3] incubation;            // [1] = T/F inclusion of time from O to A
                                             // (incubation period); [2] = mean [3] = sd
 real<lower = 0> upper_bound;                // the latest time of observation
 // [1] = mean [2] = sd
 int A;                                      // number of clusters
 array[N] int<lower=1,upper=A> cluster;
}
transformed data {
}
parameters {
 real log_mean_intercept;
 vector[A] beta_log_mean_cluster;
 real<lower=0> sigma_log_mean_cluster;
 real log_alpha_intercept;
 vector[A] beta_log_alpha_cluster;
 real<lower=0> sigma_log_alpha_cluster;
 vector<lower = 0, upper = 1>[N] a_window;   // where time a lies in the event A window
 vector<lower = 0, upper = 1>[N] b_window;   // where time b lies in the event B window
 vector<lower = 0, upper = 1>[N] a2_window;  // where time a2 lies in the event A window
 vector<lower = 0>[N] t0;                    // time from O to A
}
transformed parameters {
 vector[A] logmean_ = log_mean_intercept + beta_log_mean_cluster*sigma_log_mean_cluster; // the distribution mean
 vector[A] alpha = exp(log_alpha_intercept + beta_log_alpha_cluster*sigma_log_alpha_cluster);
 vector[A] mean_ = exp(logmean_);
 vector<lower = 0>[A] beta = mean_./tgamma(1.0 + 1.0./alpha);
 vector<lower = min(a_minus), upper = max(a_plus)>[N] a;
 vector<lower = min(a_minus), upper = max(a_plus)>[N] a2;
 vector<lower = min(b_minus), upper = max(b_plus)>[N] b;
 vector[N] ub;
 b = b_minus + (b_plus - b_minus) .* b_window;
 for (n in 1:N){
   ub[n] = min([a_plus[n], b[n]]');
 }
 a = a_minus + (ub - a_minus) .* a_window;
 a2 = a_minus + (ub - a_minus) .* a2_window;
 for (n in 1:N){
   ub[n] = min([upper_bound,a[n]+11]');
 }
}
model {
 log_mean_intercept ~ normal(1, 2);
 beta_log_mean_cluster ~ normal(0,1);
 sigma_log_mean_cluster ~ exponential(1);
 log_alpha_intercept ~ normal(1.6, 1);
 beta_log_alpha_cluster ~ normal(0,1);
 sigma_log_alpha_cluster ~ exponential(1);
 if (incubation[1] == 1){
   t0 ~ lognormal(incubation[2], incubation[3]);
 } else {
   t0 ~ normal(0, 1e-10);
 }
   target += weibull_lpdf((b - a + t0) | alpha[cluster], beta[cluster])
   - log(weibull_cdf((upper_bound - a + t0) | alpha[cluster], beta[cluster]) - weibull_cdf(1 | alpha[cluster], beta[cluster]));

}
generated quantities {
 vector[A] sd_ = beta.*sqrt(tgamma(1.0+2.0./alpha)-(tgamma(1.0+1.0./alpha)).^2);
 vector[A] limit_val_ = beta.*(-log(1-0.95)).^(1./alpha);
 vector[A] median_ = beta.*(-log(1-0.5)).^(1./alpha);
 vector[N] log_likelihood;
 for (n in 1:N) {
   log_likelihood[n] = weibull_lpdf(b[n] - a[n] + t0[n] | alpha[cluster[n]], beta[cluster[n]])
   - log(weibull_cdf((upper_bound - a[n] + t0[n]) | alpha[cluster[n]], beta[cluster[n]]) - weibull_cdf(1 | alpha[cluster[n]], beta[cluster[n]]));
 }
}
