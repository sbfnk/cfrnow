
<!-- README.md is generated from README.Rmd. Please edit that file -->

# cfrnow

<!-- badges: start -->

[![R-CMD-check](https://github.com/sbfnk/cfrnow/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/sbfnk/cfrnow/actions/workflows/R-CMD-check.yaml)
[![codecov](https://codecov.io/gh/sbfnk/cfrnow/branch/main/graph/badge.svg)](https://app.codecov.io/gh/sbfnk/cfrnow)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

Real-time case fatality ratio (CFR) estimation from line-list data,
using a Bayesian mixture-cure survival model. It is the R counterpart of
the CFR component of the Julia
[`bdbv-2026-linelist-analysis`](https://github.com/sbfnk/bdbv-2026-linelist-analysis)
model, built for collaborators who work in R.

## Why

The naive `deaths / cases` ratio underestimates in real time: recent
cases have not yet had time to die. Restricting to *resolved* cases,
`deaths / (deaths + recoveries)`, overestimates instead, because deaths
resolve faster than recoveries, so the resolved set is enriched for
deaths.

`cfrnow` conditions on neither. Each case is fatal with probability
`cfr`; a fatal case dies at an interval-censored onset-to-death delay
`F`; a case still unresolved at the cut-off is right-censored,
contributing `1 - cfr * F(t)` (it is non-fatal, or fatal but not yet
resolved). With `F` fixed this is the Ghani/Nishiura estimator; here `F`
is co-estimated and its uncertainty propagated. Recovery dates are
optional — when recorded, a `recovery_delay` adds a two-outcome fit that
also uses recovery *timing* (`(1 - cfr) f_R(r)`), sharpening the
estimate. Onset interval-censoring and real-time right-truncation of `F`
are handled by
[primarycensored](https://primarycensored.epinowcast.org/).

## Usage

``` r
library(cfrnow)
#> Loading required package: distspec
#> 
#> Attaching package: 'distspec'
#> The following objects are masked from 'package:stats':
#> 
#>     Gamma, sd

# your own onset_date / death_date line list, or a simulated one
set.seed(1)
ll <- simulate_linelist(n = 400, cfr = 0.55, onset_days = 45,
                        delay = Gamma(mean = 12.75, sd = 7))
d <- prepare_cfr_data(ll, obs_time = as.Date("2026-02-20"))
c(cases = d$n_cases, deaths = d$n_deaths, censored = d$n_cens)
#>    cases   deaths censored 
#>      400      180      220
```

Supply the onset-to-death `delay` and a `cfr_prior` as
[distspec](https://epiforecasts.io/distspec/) distributions — neither
has a default. In the delay, a native parameter can be a `Normal()`
prior (co-estimated) or a fixed number (held fixed; fixing the whole
delay gives the Ghani/Nishiura estimator). The `cfr_prior` is a `Beta()`
and matters because the CFR is weakly identified early on: `Beta(1, 1)`
is uniform, `Beta(1, 9)` favours a low CFR, and `Beta(6.6, 13.4)` (mean
0.33) suits a high-fatality pathogen. Here we use the BDBV/Isiro
onset-to-death prior:

``` r
library(distspec)
onset_to_death <- LogNormal(meanlog = Normal(2.41, 0.2), sdlog = Normal(0.51, 0.15))
fit <- fit_cfr(d, delay = onset_to_death, cfr_prior = Beta(1, 1))
summary(fit)
#>     quantity   mean   q2.5    q50  q97.5  rhat ess_bulk
#> 1        cfr  0.560  0.505  0.560  0.615 1.001     3100
#> 2 delay_mean 12.700 11.600 12.700 13.900 1.000     3300
#> 3   delay_sd  7.100  6.200  7.000  8.200 1.000     3000
```

Pass `obs_time = NULL` for a retrospective fit (it reduces to the naive
ratio with the delay as a nuisance). Swap the family or add recovery
timing:

``` r
fit_cfr(d, delay = Gamma(shape = Normal(3, 1), rate = Normal(0.25, 0.1)),
        cfr_prior = Beta(1, 1))                                           # gamma
fit_cfr(d, delay = LogNormal(meanlog = 2.41, sdlog = 0.51),
        cfr_prior = Beta(1, 1))                                           # fixed-F
fit_cfr(d, delay = onset_to_death, cfr_prior = Beta(1, 1),
        recovery_delay = LogNormal(Normal(2.9, 0.3), Normal(0.5, 0.2)))   # two-outcome
```

To estimate the delay from a line list instead, use
[epidist](https://epidist.epinowcast.org/) (which handles the double
interval censoring and real-time right-truncation) and feed it in as the
`delay` — see `vignette("cfrnow")`, which also benchmarks cfrnow against
Aalen–Johansen and Fine–Gray.

Your line list needs `onset_date` and `death_date` (`NA` for cases that
have not died). Optional: `onset_lower`/`onset_upper` widen the onset
censoring, and `recovery_date` marks a non-fatal recovery. `summary()`
reports `rhat`/`ess_bulk` and flags `cfr_low_information` when the CFR
posterior has barely moved from its prior — expected early on, when the
estimate is prior-driven.

## Assumptions and caveats

The correction fixes a *timing* bias; it cannot repair the data. Read
these before quoting a number:

- **Complete death ascertainment.** A death that never reaches the line
  list (say a community death outside a treatment centre) is treated as
  a survivor and biases the CFR down — the biggest threat where
  ascertainment is ETC-centred.
- **Use the death notification date** in `death_date`, not the true date
  of death: a case counts as a death once its `death_date` is on or
  before the cut-off. If notification lags, this lets the censoring
  absorb the delay and keeps `cfr` correct — but the delay the model
  then works with is onset-to-*notification*, so read
  `delay_mean`/`delay_sd` (and any biological delay prior) with that in
  mind. With no notification lag, plain death dates are fine. A lag in
  notifying an *onset* does not bias `cfr` unless notification depends
  on outcome.
- **Stationary delay and homogeneous CFR.** One delay and CFR over the
  whole outbreak; as treatment scales up the true CFR should fall, and a
  pooled estimate lags reality. Time-varying CFR is on the roadmap.
- **Delay family is chosen, not estimated.** Gamma and lognormal differ
  in tail weight (how fast a recent case counts as “probably cured”) and
  the family is not identifiable from sparse data, so check sensitivity
  by refitting the other one.
- **Recovery timing leans on complete discharge data.** In the
  two-outcome fit, an actually recovered case whose recovery is
  unrecorded stays censored and gets pushed toward the fatal branch over
  time, biasing `cfr` *up* — the mirror of the death-undercount bias.
  The death-only default is insensitive to this (a non-fatal case
  contributes `1 - cfr` regardless), so prefer it where discharge
  recording is patchy.

## Roadmap

Known gaps, roughly in priority order: time-varying `cfr` (e.g. a random
walk on `logit(cfr)`); a death-reporting-delay nowcasting layer;
stratification (age, sex, vaccination); and posterior-predictive checks.

## Requirements

- [cmdstanr](https://mc-stan.org/cmdstanr/) and a working CmdStan
  install
- [primarycensored](https://primarycensored.epinowcast.org/)
- [distspec](https://epiforecasts.io/distspec/)

The primarycensored Stan functions are vendored into
`inst/stan/include/pcd_functions.stan`; regenerate with
`data-raw/vendor_stan.R` after upgrading primarycensored.
