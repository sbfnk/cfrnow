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
  for (fam in curecfr_families()) {
    ll <- simulate_linelist(n = 20000, cfr = 1, delay_mean = 12.75,
                            delay_sd = 7, delay_family = fam)
    delay <- as.numeric(ll$death_date - ll$onset_date)
    expect_equal(mean(delay), 12.75, tolerance = 0.3)
    expect_equal(sd(delay), 7, tolerance = 0.3)
  }
})

test_that("unsupported families are rejected", {
  expect_error(delay_to_native(10, 5, "weibull"))
  expect_error(delay_dist_id("weibull"))
})
