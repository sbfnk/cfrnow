test_that("gamma mean/sd converts to shape/rate with the right moments", {
  m <- 12.75; s <- 7
  p <- delay_to_native(m, s, "gamma")
  # gamma(shape, rate): mean = shape / rate, sd = sqrt(shape) / rate
  expect_equal(p[["shape"]] / p[["rate"]], m)
  expect_equal(sqrt(p[["shape"]]) / p[["rate"]], s)
})

test_that("lognormal mean/sd converts to meanlog/sdlog with the right moments", {
  m <- 12.75; s <- 7
  p <- delay_to_native(m, s, "lognormal")
  expect_equal(exp(p[["meanlog"]] + p[["sdlog"]]^2 / 2), m)
  expect_equal(sqrt((exp(p[["sdlog"]]^2) - 1)) * m, s)
})

test_that("simulated delays match the requested mean and sd", {
  set.seed(1)
  for (fam in cfrnow_families()) {
    ll <- simulate_linelist(n = 20000, cfr = 1, delay_mean = 12.75,
                            delay_sd = 7, delay_family = fam)
    delay <- as.numeric(ll$death_date - ll$onset_date)
    expect_equal(mean(delay), 12.75, tolerance = 0.3)
    expect_equal(sd(delay), 7, tolerance = 0.3)
  }
})

test_that("recovery = TRUE adds recoveries for non-fatal cases only", {
  set.seed(1)
  ll <- simulate_linelist(n = 300, cfr = 0.5, recovery = TRUE)
  expect_true("recovery_date" %in% names(ll))
  fatal <- !is.na(ll$death_date)
  expect_true(all(is.na(ll$recovery_date[fatal])))     # fatal: no recovery
  expect_true(all(!is.na(ll$recovery_date[!fatal])))   # non-fatal: recovered
  expect_false("recovery_date" %in% names(simulate_linelist(n = 5)))  # off by default
})

test_that("unsupported families are rejected", {
  expect_error(delay_to_native(10, 5, "weibull"))
  expect_error(delay_native_order("weibull"))
  expect_error(stan_delay_fields(dist.spec::Normal(mean = 5, sd = 1),
                                 "dist_id", "p"))
})

test_that("stan_delay_fields reads estimated (Normal) native parameters", {
  d <- dist.spec::LogNormal(meanlog = dist.spec::Normal(2.4, 0.2),
                            sdlog = dist.spec::Normal(0.5, 0.15))
  ds <- stan_delay_fields(d, "dist_id", "p")
  expect_equal(ds$dist_id, primarycensored::pcd_stan_dist_id("lognormal", "delay"))
  expect_equal(ds$p1_est, 1L)
  expect_equal(ds$p2_est, 1L)
  expect_equal(ds$p1_prior_mean, 2.4)
  expect_equal(ds$p2_prior_sd, 0.15)
})

test_that("stan_delay_fields reads a fixed delay (fixed-F), with a q prefix too", {
  d <- dist.spec::Gamma(shape = 3, rate = 0.24)
  ds <- stan_delay_fields(d, "recovery_dist_id", "q")
  expect_equal(ds$recovery_dist_id,
               primarycensored::pcd_stan_dist_id("gamma", "delay"))
  expect_equal(ds$q1_est, 0L)
  expect_equal(ds$q2_est, 0L)
  expect_equal(ds$q1_fixed, 3)
  expect_equal(ds$q2_fixed, 0.24)
})
