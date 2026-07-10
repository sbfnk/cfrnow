# cfrnow 0.1.0

First release.

* `fit_cfr()` estimates a real-time case fatality ratio from line-list data with a
  Bayesian mixture-cure survival model. It is registered as an `epidist` model
  type, so the CFR and the onset-to-death delay both take `brms` formulas.
* `prepare_cfr_data()` turns a line list into model inputs. It sorts each case, at
  a chosen observation cut-off, into an observed death, a resolved non-death, or a
  right-censored survivor.
* The onset-to-death delay (LogNormal or Gamma) can be co-estimated or held fixed.
  Hold it fixed and you get the Ghani/Nishiura estimator.
* Pass a `recovery_date` column and a two-outcome fit also times recoveries.
* Put a `brms` formula on the CFR or the delay for covariates or a time-varying
  CFR.
* `simulate_linelist()` builds line lists for testing and examples.
* `summary()` and `print()` report the corrected CFR, the delay moments,
  convergence diagnostics, and a flag for when the CFR is only weakly identified.
