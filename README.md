
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

The naive `deaths / cases` ratio is biased downward in real time: recent
cases have not yet had time to die. Restricting to cases with a
*resolved* outcome (`deaths / (deaths + recoveries)`) swaps that for an
upward bias, because deaths resolve faster than recoveries, so the
resolved set is enriched for deaths at any mid-outbreak cut-off.

`cfrnow` avoids conditioning on resolution at all. Each case is fatal
with probability `cfr`; a fatal case dies at an interval-censored
onset-to-death delay `F`; every case still alive at the cut-off is
right-censored, contributing the mixture-cure survival term
`1 - cfr * F(t)`, the probability that it is either non-fatal or fatal
but not yet resolved. With `F` fixed this is the Ghani/Nishiura
`deaths / sum_i F(t_i)` estimator; here `F` is co-estimated and its
uncertainty propagated. The enrichment bias is avoided because
unresolved cases are right-censored, not dropped. So cfrnow needs only
death timing and who is still a case, and does not depend on recovery
dates being recorded. When they are recorded, supplying a
`recovery_delay` fits a competing-risks model that also uses the
recovery *timing* (`(1 - cfr) f_R(r)` for a recovery at `r`), which
sharpens the estimate.

Onset dates are interval-censored and, in real-time mode, `F` is
right-truncated at the cut-off. Both are handled by the analytical
censored-CDF machinery vendored from
[`primarycensored`](https://primarycensored.epinowcast.org/).

## Usage

``` r
library(cfrnow)

# your own onset_date / death_date line list, or a simulated one
ll <- simulate_linelist(n = 400, cfr = 0.55, onset_days = 45)
d <- prepare_cfr_data(ll, obs_time = as.Date("2026-02-20"))
c(cases = d$n_cases, deaths = d$n_deaths, censored = d$n_cens)
#>    cases   deaths censored 
#>      400      195      205
```

``` r
fit <- fit_cfr(d)
summarise_cfr(fit)
#>     quantity   mean   q2.5    q50  q97.5  rhat ess_bulk
#> 1        cfr  0.560  0.505  0.560  0.615 1.001     3100
#> 2 delay_mean 12.700 11.600 12.700 13.900 1.000     3300
#> 3   delay_sd  7.100  6.200  7.000  8.200 1.000     3000
```

Pass `obs_time = NULL` for a retrospective (fully-resolved) fit, which
reduces to the naive ratio with the delay as a nuisance.

The onset-to-death delay is given as a
[dist.spec](https://epiforecasts.io/dist.spec/) distribution via the
`delay` argument. `default_delay()` is a lognormal with `Normal()`
priors on its parameters; supply a `dist.spec::Gamma()` to change
family, or give a parameter as a number / `dist.spec::Fixed()` to hold
it fixed:

``` r
library(dist.spec)
fit_cfr(d, delay = Gamma(shape = Normal(3, 1), rate = Normal(0.25, 0.1)))  # gamma
fit_cfr(d, delay = LogNormal(meanlog = 2.41, sdlog = 0.51))                # fixed-F
fit_cfr(d, recovery_delay = default_recovery_delay())  # competing risks (recovery timing)
```

Which of these to use follows from your data: onset and death dates give
the default death-only fit; reliable recovery dates let you add a
`recovery_delay` for the competing-risks fit; and a delay you trust from
elsewhere can be fixed for the Ghani/Nishiura estimator.

Your line list needs an `onset_date` column and a `death_date` column
(`NA` for cases that have not died). Optional columns:
`onset_lower`/`onset_upper` give an onset window that widens the primary
censoring, and `recovery_date` marks a non-fatal recovery.

`summarise_cfr()` reports `rhat`/`ess_bulk` and flags
(`cfr_low_information`) when the CFR posterior has barely moved from its
prior. That is expected early in an outbreak, when few deaths have
resolved and the CFR is only weakly identified. Treat a flagged estimate
as prior-driven.

## Assumptions and caveats

The correction fixes a *timing* bias; it cannot repair the data. Read
these before quoting a number:

- **Complete death ascertainment.** A death that never reaches the line
  list (for example a community death outside a treatment centre) is
  silently treated as a survivor, biasing the CFR down. In an outbreak
  where ascertainment is ETC-centred this is the single biggest threat
  to validity.
- **Deaths known on their day of occurrence.** The model assumes a death
  enters the data when it happens. Real notification lag reintroduces
  the very bias the model exists to remove, through the data pipeline
  rather than the delay. If that lag is non-trivial, nowcast or caveat
  the death series.
- **Stationary delay and homogeneous CFR.** A single onset-to-death
  delay and CFR over the whole outbreak. As treatment access (supportive
  care, monoclonals) scales up, the true CFR should fall, and a pooled
  estimate lags reality precisely when that matters. Time-varying CFR is
  on the roadmap below.
- **Delay family is chosen, not estimated.** Gamma and lognormal differ
  in tail weight, which affects how quickly a recent case counts as
  “probably cured”. With sparse data the family is not identifiable from
  the line list, so check sensitivity by refitting with a
  `dist.spec::Gamma()` delay.
- **Recovery timing leans on discharge data.** The competing-risks fit
  is only as good as the recovery dates and the recovery-delay prior;
  where discharge recording is patchy, the death-only default is the
  safer choice.

## Roadmap

Known gaps, roughly in priority order: time-varying `cfr` (e.g. a random
walk on `logit(cfr)`); a death-reporting-delay nowcasting layer;
stratification (age, sex, vaccination); and posterior-predictive checks.

## Status

Early but functional. The model, data preparation, summaries, a vignette
and unit tests are in place, and `R CMD check` is clean. Lognormal (the
default) and gamma onset-to-death delays are supported via
[dist.spec](https://epiforecasts.io/dist.spec/), with recovery timing
available through a competing-risks fit.

## Requirements

- [cmdstanr](https://mc-stan.org/cmdstanr/) and a working CmdStan
  install
- [primarycensored](https://primarycensored.epinowcast.org/)
- [dist.spec](https://epiforecasts.io/dist.spec/)

The primarycensored Stan functions are vendored into
`inst/stan/include/pcd_functions.stan`; regenerate with
`data-raw/vendor_stan.R` after upgrading primarycensored.
