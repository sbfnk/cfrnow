# Coverage without a sampler: drive model generation through make_stancode,
# which exercises the family / formula / stancode / translate paths (cure_model.R
# and translate.R) but needs no CmdStan, so it runs on CI where the fit tests skip.

stancode_for <- function(cure, delay, recovery_delay = NULL) {
  dd <- .delay_family_prior(delay)
  prior <- c(.cfr_prior_to_brms(Beta(1, 1)), dd$prior)
  if (!is.null(recovery_delay) && isTRUE(attr(cure, "use_recovery"))) {
    rd <- .delay_family_prior(recovery_delay, main = FALSE)
    prior <- c(prior, rd$prior)
    attr(cure, "recovery_family") <- brms:::validate_family(rd$family)
  }
  epidist::epidist(cure,
    formula = mu ~ 1, family = dd$family, prior = prior,
    merge_priors = FALSE, fn = brms::make_stancode
  )
}

test_that("a death-only lognormal fit generates the cure lpmf", {
  cure <- as_epidist_cure_model(prepare_cfr_data(
    simulate_linelist(n = 100, cfr = 0.4, delay = LogNormal(2.4, 0.5)),
    obs_time = NULL
  ))
  code <- stancode_for(cure, LogNormal(meanlog = 2.41, sdlog = 0.51))
  expect_true(grepl("cfrnow_lognormal_lpmf", code))
  expect_true(grepl("primarycensored", code))
})

test_that("a gamma fit generates a gamma cure lpmf", {
  cure <- as_epidist_cure_model(prepare_cfr_data(
    simulate_linelist(n = 100, cfr = 0.4, delay = Gamma(mean = 8, sd = 4)),
    obs_time = NULL
  ))
  code <- stancode_for(cure, Gamma(shape = Normal(4, 1), rate = Normal(0.5, 0.2)))
  expect_true(grepl("cfrnow_gamma_lpmf", code))
})

test_that("a two-outcome fit generates recovery branches with its own family", {
  ll <- simulate_linelist(
    n = 200, cfr = 0.4, delay = Gamma(mean = 12, sd = 6),
    recovery = LogNormal(mean = 20, sd = 8)
  )
  cure <- as_epidist_cure_model(prepare_cfr_data(ll, obs_time = NULL))
  expect_true(attr(cure, "use_recovery"))
  code <- stancode_for(
    cure, Gamma(shape = Normal(4, 1), rate = Normal(0.3, 0.1)),
    recovery_delay = LogNormal(meanlog = Normal(2.9, 0.3), sdlog = Normal(0.5, 0.2))
  )
  expect_true(grepl("outcome == 2", code)) # timed-recovery branch present
  expect_true(grepl("rmu", code)) # recovery params are r-prefixed
})

test_that("an intercept-free cfr formula routes the prior to the coefficients", {
  # where cfr_prior lands depends on whether the cfr formula keeps its intercept
  expect_true(.cfr_has_intercept(mu ~ 1)) # cfr defaults to intercept-only
  expect_true(.cfr_has_intercept(brms::bf(mu ~ 1, cfr ~ grp)))
  expect_false(.cfr_has_intercept(brms::bf(mu ~ 1, cfr ~ 0 + grp)))

  expect_equal(.cfr_prior_to_brms(Beta(1, 1))$class, "Intercept")
  expect_equal(.cfr_prior_to_brms(Beta(1, 1), "b")$class, "b")

  a <- simulate_linelist(n = 80, cfr = 0.3, delay = LogNormal(2.4, 0.5))
  b <- simulate_linelist(n = 80, cfr = 0.6, delay = LogNormal(2.4, 0.5))
  ca <- as_epidist_cure_model(prepare_cfr_data(a, obs_time = NULL))
  cb <- as_epidist_cure_model(prepare_cfr_data(b, obs_time = NULL))
  ca$grp <- "x"
  cb$grp <- "y"
  cure <- as_epidist_cure_model(rbind(ca, cb))

  dd <- .delay_family_prior(LogNormal(meanlog = 2.41, sdlog = 0.51))
  f <- brms::bf(mu ~ 1, cfr ~ 0 + grp)

  # the prior on the (absent) intercept is what brms rejects ...
  expect_error(
    epidist::epidist(cure,
      formula = f, family = dd$family,
      prior = c(.cfr_prior_to_brms(Beta(1, 1), "Intercept"), dd$prior),
      merge_priors = FALSE, fn = brms::make_stancode
    )
  )
  # ... and moving it onto the coefficients generates cleanly
  code <- epidist::epidist(cure,
    formula = f, family = dd$family,
    prior = c(.cfr_prior_to_brms(Beta(1, 1), "b"), dd$prior),
    merge_priors = FALSE, fn = brms::make_stancode
  )
  expect_true(grepl("cfr", code))
})
