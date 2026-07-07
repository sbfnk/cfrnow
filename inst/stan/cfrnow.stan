// Real-time case fatality ratio via a Bayesian mixture-cure survival model.
//
// Each case is fatal with probability `cfr`; a fatal case dies at an
// onset-to-death delay `F_D`. A non-fatal case is "cured" and, when recovery is
// modelled (`use_recovery` = 1), recovers at an onset-to-recovery delay `F_R`.
// At the observation cut-off a case is one of:
//   - an observed death at delay d      -> cfr * f_D(d)
//   - an observed recovery at delay r   -> (1 - cfr) * f_R(r)   [use_recovery]
//   - a resolved non-death, time unknown -> (1 - cfr)           [retrospective]
//   - still unresolved at elapsed time t ->
//       use_recovery: cfr * (1 - F_D(t)) + (1 - cfr) * (1 - F_R(t))   (two-outcome mixture)
//       otherwise:    1 - cfr * F_D(t)                                (death only)
// With the delays fixed and recovery off this is the Ghani/Nishiura estimator.
//
// Each delay is carried in its native parameters (dist.spec-style): lognormal
// (meanlog, sdlog) or gamma (shape, rate). Each native parameter is either
// fixed (supplied as data) or estimated (a Normal prior); the length-0/1 arrays
// declare a sampled parameter only when it is estimated, so one model covers
// fixed, estimated and mixed delays. Onset is interval-censored (a per-case
// day-window via primarycensored's primary uniform); deaths and recoveries are
// censored to their recorded day; the cut-off is an exact boundary, so
// unresolved cases are only primary-censored.
functions {
#include include/pcd_functions.stan

  // Mean and sd (days) of a native (p1, p2) pair, family chosen by dist_id.
  vector delay_moments(int did, real a, real b) {
    vector[2] m;
    if (did == 1) {              // lognormal: (meanlog, sdlog)
      m[1] = exp(a + square(b) / 2);
      m[2] = m[1] * sqrt(expm1(square(b)));
    } else {                     // gamma: (shape, rate)
      m[1] = a / b;
      m[2] = sqrt(a) / b;
    }
    return m;
  }
}
data {
  // Deaths observed by the cut-off: integer onset->death delay + onset window.
  int<lower=0> n_death;
  array[n_death] int<lower=0> death_delay;
  vector<lower=0>[n_death] death_width;
  // Recoveries observed by the cut-off: integer onset->recovery delay + window.
  int<lower=0> n_recovery;
  array[n_recovery] int<lower=0> recovery_delay;
  vector<lower=0>[n_recovery] recovery_width;
  // Non-deaths resolved but with no recovery time (retrospective): log(1-cfr).
  int<lower=0> n_resolved;
  // Cases still unresolved at a finite cut-off: elapsed follow-up + onset window.
  int<lower=0> n_cens;
  vector<lower=0>[n_cens] censor_time;
  vector<lower=0>[n_cens] censor_width;
  // Whether to model the onset-to-recovery delay (two-outcome mixture).
  int<lower=0, upper=1> use_recovery;
  int<lower=1> primary_id;
  // Onset-to-death delay family + native parameters (p1, p2).
  int<lower=1> dist_id;
  int<lower=0, upper=1> p1_est;
  int<lower=0, upper=1> p2_est;
  real<lower=0> p1_fixed;
  real<lower=0> p2_fixed;
  real p1_prior_mean;
  real<lower=0> p1_prior_sd;
  real p2_prior_mean;
  real<lower=0> p2_prior_sd;
  // Onset-to-recovery delay family + native parameters (q1, q2). Ignored when
  // use_recovery = 0 (pass any valid dummy).
  int<lower=1> recovery_dist_id;
  int<lower=0, upper=1> q1_est;
  int<lower=0, upper=1> q2_est;
  real<lower=0> q1_fixed;
  real<lower=0> q2_fixed;
  real q1_prior_mean;
  real<lower=0> q1_prior_sd;
  real q2_prior_mean;
  real<lower=0> q2_prior_sd;
  // CFR prior.
  real<lower=0> cfr_a;
  real<lower=0> cfr_b;
}
parameters {
  array[p1_est] real<lower=0> p1_par;
  array[p2_est] real<lower=0> p2_par;
  array[use_recovery * q1_est] real<lower=0> q1_par;
  array[use_recovery * q2_est] real<lower=0> q2_par;
  real<lower=0, upper=1> cfr;
}
transformed parameters {
  real p1 = p1_est == 1 ? p1_par[1] : p1_fixed;
  real p2 = p2_est == 1 ? p2_par[1] : p2_fixed;
  array[2] real params = {p1, p2};
  real q1 = (use_recovery == 1 && q1_est == 1) ? q1_par[1] : q1_fixed;
  real q2 = (use_recovery == 1 && q2_est == 1) ? q2_par[1] : q2_fixed;
  array[2] real rparams = {q1, q2};
}
model {
  array[0] real primary_params = rep_array(0.0, 0);  // uniform primary: no params

  if (p1_est == 1) {
    p1_par[1] ~ normal(p1_prior_mean, p1_prior_sd);
  }
  if (p2_est == 1) {
    p2_par[1] ~ normal(p2_prior_mean, p2_prior_sd);
  }
  if (use_recovery == 1 && q1_est == 1) {
    q1_par[1] ~ normal(q1_prior_mean, q1_prior_sd);
  }
  if (use_recovery == 1 && q2_est == 1) {
    q2_par[1] ~ normal(q2_prior_mean, q2_prior_sd);
  }
  cfr ~ beta(cfr_a, cfr_b);

  // Every non-death (untimed resolved + observed recovery) carries the cure
  // factor (1 - cfr); observed recoveries additionally carry the recovery-delay
  // density when recovery is modelled.
  target += (n_resolved + n_recovery) * log1m(cfr);
  if (use_recovery == 1) {
    for (j in 1:n_recovery) {
      target += primarycensored_lpmf(recovery_delay[j] | recovery_dist_id,
          rparams, recovery_width[j], recovery_delay[j] + 1.0,
          0.0, positive_infinity(), primary_id, primary_params);
    }
  }

  // Observed deaths: each fatal (log cfr) with a doubly-interval-censored delay.
  target += n_death * log(cfr);
  for (i in 1:n_death) {
    target += primarycensored_lpmf(death_delay[i] | dist_id, params,
        death_width[i], death_delay[i] + 1.0,
        0.0, positive_infinity(), primary_id, primary_params);
  }

  // Unresolved cases at the cut-off.
  for (i in 1:n_cens) {
    real fbar_d = primarycensored_cdf(censor_time[i] | dist_id, params,
        censor_width[i], 0.0, positive_infinity(), primary_id, primary_params);
    if (use_recovery == 1) {
      // Two-outcome mixture: P(fatal not yet dead) + P(non-fatal not yet recovered).
      real fbar_r = primarycensored_cdf(censor_time[i] | recovery_dist_id,
          rparams, censor_width[i], 0.0, positive_infinity(),
          primary_id, primary_params);
      target += log_sum_exp(log(cfr) + log1m(fbar_d),
                            log1m(cfr) + log1m(fbar_r));
    } else {
      // Death only: non-fatal cases assumed not to resolve.
      target += log1m(cfr * fbar_d);
    }
  }
}
generated quantities {
  vector[2] md = delay_moments(dist_id, p1, p2);
  vector[2] mr = delay_moments(recovery_dist_id, q1, q2);
  real delay_mean = md[1];
  real delay_sd = md[2];
  real recovery_mean = mr[1];
  real recovery_sd = mr[2];
}
