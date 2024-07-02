functions {
   // Define a function to compute the log PDF of the Bates distribution
  real bates_log_pdf(real x, int n) {
    // Compute the log PDF as the sum of the log PDFs of n standard normal distributions
    real log_pdf = 0;
    for (i in 1:n) {
      log_pdf = log_pdf + uniform_lpdf(x / sqrt(n)|0, 1); // Assuming standard normal distributions
    }
    // Add the normalisation term in log space
    log_pdf = log_pdf - 0.5 * log(2 * pi() * n);
    return log_pdf;
  }



   
    vector sens_given_pos_func(vector dfso, real shift, vector rate, real ub, vector lb, int N) {

      vector[N] y;

      y = inv_logit(shift - dfso .* rate) .* (ub -lb) + lb;

      return y;
    }

    // vector sens_given_pos_func_alt(vector dfso, real shift, real rate, real ub, real lb, int N) {
    vector sens_given_pos_func_alt(vector dfso, real shift, real rate, real ub, real lb, int N) {

      vector[N] y;

      y = inv_logit(shift - dfso .* rate) .* (ub -lb) + lb;

      return y;
    }


}


data {

  int<lower=1> A; //Num age groups

  //////////////////////////////////////
  // Duration of Positivity data
  int<lower = 1> N_dur_pos;                 // number unique (A minus, A minus, B minus, B plus) combos
  array[N_dur_pos] int<lower = 1> count;              // how many times the nth event combo is observed
  vector<lower = 0>[N_dur_pos] a_minus;            // lower limit of event A
  vector<lower = 0>[N_dur_pos] a_plus;             // upper limit of event A
  vector<lower = 0>[N_dur_pos] b_minus;            // lower limit of event B
  vector<lower = 0>[N_dur_pos] b_plus;             // upper limit of event B
  vector<lower = 0, upper = 3>[N_dur_pos] censor;  // censoring parameter (0 for no censoring, 1 for right-censoring, 2 for left-censoring, 3 for left-right censored)
  array[N_dur_pos] int<lower=1, upper=A> age_group_int;
  int maxT_dur;

  //////////////////////////////////////
  // sensitivity given positive (sgp) data
  int<lower=0> N_sgp_obs;
  array[N_sgp_obs] int N_tests;
  array[N_sgp_obs] int N_pos;
  array[N_sgp_obs] int sgp_age_group_int;
  vector[N_sgp_obs] sgp_dfso;

  // used when calculating the sens given pos functions
  int<lower=0> N_dfso;
  vector[N_dfso] dfso;

}

transformed data {

  vector[maxT_dur] dfso_conv;
  for(n in 1:maxT_dur){
    dfso_conv[n] = n-1;
  }


}

parameters {

  //////////////////////////////////////////////////////////////////////////////
  // Duration of Positivity parameters
  real dur_pos_mu_intercept;

  vector[A] dur_pos_mu_age;
  real<lower=0> dur_pos_sigma;
  real<lower=0> dur_pos_sigma_age;

  vector<lower = 0, upper = 1>[N_dur_pos] a_window; // where time a lies in the event A window
  vector<lower = 0, upper = 1>[N_dur_pos] b_window; // where time b lies in the event B window

  //////////////////////////////////////////////////////////////////////////////
  // Sensitivity given positivity parameters
  real<lower=0, upper=1> lb_mu;
  real<lower=0> lb_precision_raw;
  real<lower=0> rate_mu;
  vector[A] rate_beta;
  real<lower=0> rate_sigma;
 
  real shift_raw;
  vector<lower=0, upper=1>[A] lower_bound;
  real<lower=0, upper=1> upper_bound;
}

transformed parameters {

  ///////////////////////////////////
  // Duration of Positivity
  vector<lower = min(a_minus), upper = max(a_plus)>[N_dur_pos] a;
  vector<lower = min(b_minus), upper = max(b_plus)>[N_dur_pos] b;
  vector[N_dur_pos] ub;
  real<lower = 0> mean_val;
  real<lower = 0> sd_val;

  vector[A] dur_pos_mu;

  mean_val = exp((2.5 * dur_pos_mu_intercept + 0.5) + (dur_pos_sigma^2)/2);
  sd_val = mean_val * sqrt(exp(dur_pos_sigma^2) - 1);

  for (i in 1:A){
    dur_pos_mu[i] = (2.5 * dur_pos_mu_intercept + 0.5) + dur_pos_mu_age[i] * dur_pos_sigma_age;
  }

  // duration of positivity bounds calculation
  b = b_minus + (b_plus - b_minus) .* b_window;
  for (n in 1:N_dur_pos) {
    ub[n] = min([a_plus[n], b[n]]');
  }
  a = a_minus + (ub - a_minus) .* a_window;

  // Duration Simplex
  array[A] vector<lower=0,upper=1>[maxT_dur] prob_still_pos;
  for (i in 1:A){
    for(t in 1:maxT_dur){
      prob_still_pos[i,t] = 1-lognormal_cdf(t-1 | dur_pos_mu[i], dur_pos_sigma);
    }
  }

  // Duration pdf Simplex
  array[A] vector<lower=0>[maxT_dur*10] pos_dur_pdf;
  for (i in 1:A){
    for(t in 1:maxT_dur*10){
      pos_dur_pdf[i,t] = exp(lognormal_lpdf((t-1)/10.0 | dur_pos_mu[i], dur_pos_sigma));
    }
  }


  //////////////////////////////////////////////////////////////////////////////
  // sens given positivity

  real shift = 4 + sqrt(3) * shift_raw;
  real<lower=0> lb_precision = 100 * lb_precision_raw;
  real lb_alpha = lb_mu * lb_precision;
  real lb_beta = (1-lb_mu) * lb_precision;

  vector[A] rate = 1/(rate_mu + rate_beta * rate_sigma);


  vector<lower=0, upper=1>[N_sgp_obs] sens_given_pos_fitting = sens_given_pos_func(
    sgp_dfso,
    shift,
    rate[sgp_age_group_int],
    upper_bound,
    lower_bound[sgp_age_group_int],
    N_sgp_obs);

  array[A] vector[maxT_dur] sens_given_pos_conv;

  for(age_group_conv in 1:A){
    sens_given_pos_conv[age_group_conv] = sens_given_pos_func_alt(
      dfso,
      shift,
      rate[age_group_conv],
      upper_bound,
      lower_bound[age_group_conv],
      maxT_dur
    );

  }
 
  array[A] vector<lower=0, upper=1>[maxT_dur] sens_curve;

  for (i in 1:A){
   sens_curve[i] = sens_given_pos_conv[i] .* prob_still_pos[i];
  }


}

model {

  //////////////////////////////////////////////////////////////////////////////
  // Duration of positivity model
  dur_pos_mu_intercept ~ std_normal();  // implicit N(0.5, 0.5)
  dur_pos_mu_age ~ normal(0,1);
  dur_pos_sigma ~ inv_gamma(10, 2);
  dur_pos_sigma_age ~ exponential(5);

  for (i in 1:N_dur_pos){
    target += count[i]*beta_lpdf(a_window[i]|(12*count[i]-4)/8.0,(12*count[i]-4)/8.0);
    target += count[i]*beta_lpdf(b_window[i]|(12*count[i]-4)/8.0,(12*count[i]-4)/8.0);
  }

  for (i in 1:N_dur_pos) {

    if(censor[i] == 0){
      // no censoring
      target += count[i] * lognormal_lpdf(b[i] - a[i] | dur_pos_mu[age_group_int[i]], dur_pos_sigma);

    } else if (censor[i] == 1) {
      // right censoring
      target += count[i] * lognormal_lccdf(b_minus[i] - a[i] + 0.1| dur_pos_mu[age_group_int[i]], dur_pos_sigma);

    } else if (censor[i] == 2) {
      // left censoring
      target += count[i] * lognormal_lccdf(b[i] - a_plus[i] | dur_pos_mu[age_group_int[i]], dur_pos_sigma);

    } else if (censor[i] == 3) {
      // left-right censoring
      target += count[i] * lognormal_lccdf(b_minus[i] - a_plus[i] | dur_pos_mu[age_group_int[i]], dur_pos_sigma);
    }

  }

  // mean_val_age ~ normal(10, 2);

  // //////////////////////////////////////////////////////////////////////////////
  // // Sensitivity given positive model
  //
  shift_raw ~ std_normal(); // implicit N(4,3)

  rate_mu ~ normal(1.5, 1);
  rate_beta ~ normal(0,1);
  rate_sigma ~ exponential(1);

  lb_mu ~ beta(4,6);
  lb_precision_raw ~ exponential(1); // implicit exp(0.05)
  lower_bound ~ beta(lb_alpha, lb_beta);

  upper_bound ~ beta(15, 1);



  N_pos ~ binomial(N_tests, sens_given_pos_fitting);


}

generated quantities {

}
