// Real-time case fatality ratio via a Bayesian mixture-cure survival model.
//
// Each case is fatal with probability `cfr`; a fatal case dies at an
// onset-to-death delay drawn from the parametric family `dist_id`. Cases still
// alive at a finite observation cut-off are right-censored (either non-fatal,
// or fatal but not yet resolved), which corrects the naive deaths / cases
// ratio in real time. With the delay held fixed this reduces to the
// Ghani/Nishiura estimator cfr = deaths / sum_i F(t_i); here the delay is
// co-estimated and its uncertainty propagated.
//
// The delay is parameterised by its mean and standard deviation (both
// interpretable in days and family-independent) and converted to each family's
// native parameters below, so priors and posterior summaries do not depend on
// the family. Onset is interval-censored (a per-case day-window handled by
// primarycensored's primary uniform); death is censored to its recorded day;
// the cut-off is an exact boundary, so survivors are only primary-censored.
functions {
#include include/pcd_functions.stan

  // Convert delay (mean, sd) to the native two-parameter vector expected by
  // primarycensored for the given family. Extend the branch to add families.
  array[] real delay_to_native(real m, real s, int dist_id) {
    array[2] real p;
    if (dist_id == 1) {             // lognormal: (meanlog, sdlog)
      real cv2 = square(s / m);
      p[1] = log(m) - 0.5 * log1p(cv2);
      p[2] = sqrt(log1p(cv2));
    } else if (dist_id == 2) {      // gamma: (shape, rate)
      p[1] = square(m / s);
      p[2] = m / square(s);
    } else {
      reject("unsupported dist_id for the mean/sd parameterisation: ", dist_id);
    }
    return p;
  }
}
data {
  // Deaths observed by the cut-off: integer onset->death delay + onset window.
  int<lower=0> n_death;
  array[n_death] int<lower=0> death_delay;
  vector<lower=0>[n_death] death_width;
  // Survivors alive at a finite cut-off: elapsed follow-up + onset window.
  int<lower=0> n_cens;
  vector<lower=0>[n_cens] censor_time;
  vector<lower=0>[n_cens] censor_width;
  // Fully-resolved non-fatal cases (retrospective fit); contribute log(1-cfr).
  int<lower=0> n_resolved;
  // Delay / primary family ids (from primarycensored::pcd_stan_dist_id()).
  int<lower=1> dist_id;
  int<lower=1> primary_id;
  // Priors (on the interpretable delay mean/sd in days, and the CFR).
  real<lower=0> delay_mean_mean;
  real<lower=0> delay_mean_sd;
  real<lower=0> delay_sd_mean;
  real<lower=0> delay_sd_sd;
  real<lower=0> cfr_a;
  real<lower=0> cfr_b;
}
parameters {
  real<lower=0> delay_mean;   // onset-to-death mean (days)
  real<lower=0> delay_sd;     // onset-to-death sd (days)
  real<lower=0, upper=1> cfr;
}
transformed parameters {
  array[2] real params = delay_to_native(delay_mean, delay_sd, dist_id);
}
model {
  array[0] real primary_params = rep_array(0.0, 0);  // uniform primary: no params

  delay_mean ~ normal(delay_mean_mean, delay_mean_sd);  // T[0,] via <lower=0>
  delay_sd ~ normal(delay_sd_mean, delay_sd_sd);        // T[0,] via <lower=0>
  cfr ~ beta(cfr_a, cfr_b);

  // Fully-resolved non-deaths: the cure component only.
  target += n_resolved * log1m(cfr);

  // Observed deaths: fatal (log cfr) x doubly-interval-censored onset->death.
  for (i in 1:n_death) {
    target += log(cfr)
      + primarycensored_lpmf(death_delay[i] | dist_id, params,
          death_width[i], death_delay[i] + 1.0,
          0.0, positive_infinity(), primary_id, primary_params);
  }

  // Right-censored survivors: mixture-cure survival 1 - cfr * Fbar(t), with
  // Fbar the onset-window-averaged (primary-censored) CDF at the cut-off.
  for (i in 1:n_cens) {
    real fbar = primarycensored_cdf(censor_time[i] | dist_id, params,
        censor_width[i], 0.0, positive_infinity(), primary_id, primary_params);
    target += log1m(cfr * fbar);
  }
}
