// Real-time case fatality ratio via a Bayesian mixture-cure survival model.
//
// Each case is fatal with probability `cfr`; a fatal case dies at an
// onset-to-death delay from the parametric family `dist_id`. Cases still alive
// at a finite observation cut-off are right-censored (either non-fatal, or
// fatal but not yet resolved), which corrects the naive deaths / cases ratio in
// real time. With the delay fixed this is the Ghani/Nishiura estimator; with
// priors on the delay it is co-estimated.
//
// The delay is carried in its native parameters (dist.spec-style): lognormal
// (meanlog, sdlog) or gamma (shape, rate). Each native parameter is either
// fixed (supplied as data) or estimated (a Normal prior); the length-0/1 arrays
// declare a sampled parameter only when it is estimated, so one model covers
// fully fixed, fully estimated, and mixed delays. `delay_mean`/`delay_sd` are
// recovered as generated quantities so summaries stay family-independent.
//
// Onset is interval-censored (a per-case day-window via primarycensored's
// primary uniform); death is censored to its recorded day; the cut-off is an
// exact boundary, so survivors are only primary-censored.
functions {
#include include/pcd_functions.stan
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
  // Native delay parameters. Each is either fixed (`*_est` = 0, value in
  // `*_fixed`) or estimated (`*_est` = 1, Normal prior in `*_prior_*`). p1/p2
  // are (meanlog, sdlog) for lognormal and (shape, rate) for gamma.
  int<lower=0, upper=1> p1_est;
  int<lower=0, upper=1> p2_est;
  real<lower=0> p1_fixed;
  real<lower=0> p2_fixed;
  real p1_prior_mean;
  real<lower=0> p1_prior_sd;
  real p2_prior_mean;
  real<lower=0> p2_prior_sd;
  // CFR prior.
  real<lower=0> cfr_a;
  real<lower=0> cfr_b;
}
parameters {
  // Sampled only when estimated. Both native parameters are positive here; for
  // lognormal this constrains meanlog > 0, harmless for onset-to-death delays
  // (a median below one day is implausible) and it keeps a single model.
  array[p1_est] real<lower=0> p1_par;
  array[p2_est] real<lower=0> p2_par;
  real<lower=0, upper=1> cfr;
}
transformed parameters {
  real p1 = p1_est == 1 ? p1_par[1] : p1_fixed;
  real p2 = p2_est == 1 ? p2_par[1] : p2_fixed;
  array[2] real params = {p1, p2};
}
model {
  array[0] real primary_params = rep_array(0.0, 0);  // uniform primary: no params

  if (p1_est == 1) {
    p1_par[1] ~ normal(p1_prior_mean, p1_prior_sd);
  }
  if (p2_est == 1) {
    p2_par[1] ~ normal(p2_prior_mean, p2_prior_sd);
  }
  cfr ~ beta(cfr_a, cfr_b);

  // Fully-resolved non-deaths: the cure component only.
  target += n_resolved * log1m(cfr);

  // Observed deaths: each fatal (log cfr) with a doubly-interval-censored delay.
  target += n_death * log(cfr);
  for (i in 1:n_death) {
    target += primarycensored_lpmf(death_delay[i] | dist_id, params,
        death_width[i], death_delay[i] + 1.0,
        0.0, positive_infinity(), primary_id, primary_params);
  }

  // Right-censored survivors: mixture-cure survival 1 - cfr * Fbar(t).
  for (i in 1:n_cens) {
    real fbar = primarycensored_cdf(censor_time[i] | dist_id, params,
        censor_width[i], 0.0, positive_infinity(), primary_id, primary_params);
    target += log1m(cfr * fbar);
  }
}
generated quantities {
  // Delay mean and sd (days): family-independent summaries of the native pair.
  real delay_mean;
  real delay_sd;
  if (dist_id == 1) {           // lognormal: p1 = meanlog, p2 = sdlog
    delay_mean = exp(p1 + square(p2) / 2);
    delay_sd = delay_mean * sqrt(expm1(square(p2)));
  } else {                      // gamma: p1 = shape, p2 = rate
    delay_mean = p1 / p2;
    delay_sd = sqrt(p1) / p2;
  }
}
