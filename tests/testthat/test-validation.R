# Input-validation guards and pure-R helpers. None of these need CmdStan, so
# they run everywhere and exercise the error paths fit_cfr() takes before it
# ever reaches the sampler.

test_that("fit_cfr rejects data that did not come from prepare_cfr_data", {
  expect_error(
    fit_cfr(data.frame(x = 1), delay = distspec::LogNormal(2.4, 0.5),
            cfr_prior = distspec::Beta(1, 1)),
    "prepare_cfr_data"
  )
})

test_that("fit_cfr requires a delay and a cfr_prior", {
  d <- prepare_cfr_data(
    simulate_linelist(n = 5, delay = distspec::LogNormal(2.4, 0.5)),
    obs_time = as.Date("2026-02-01")
  )
  expect_error(fit_cfr(d, cfr_prior = distspec::Beta(1, 1)), "delay")
  expect_error(fit_cfr(d, delay = distspec::LogNormal(2.4, 0.5)), "cfr_prior")
})

test_that("parse_delay_param handles numbers, Fixed(), Normal() and rejects rest", {
  expect_equal(parse_delay_param(3, "p")$est, 0L)
  expect_equal(parse_delay_param(3, "p")$fixed, 3)
  expect_equal(parse_delay_param(distspec::Fixed(4), "p")$fixed, 4)
  expect_equal(parse_delay_param(distspec::Fixed(4), "p")$est, 0L)
  n <- parse_delay_param(distspec::Normal(2, 0.3), "p")
  expect_equal(n$est, 1L)
  expect_equal(c(n$prior_mean, n$prior_sd), c(2, 0.3))
  expect_error(parse_delay_param(distspec::Gamma(3, 1), "p"), "Normal")
  expect_error(parse_delay_param(list(), "p"), "unrecognised")
})

test_that("naive_cfr is deaths/cases, and NA when there are no cases", {
  expect_equal(naive_cfr(list(n_cases = 10, n_deaths = 3)), 0.3)
  expect_true(is.na(naive_cfr(list(n_cases = 0, n_deaths = 0))))
})

test_that("beta_sd matches the closed-form Beta standard deviation", {
  expect_equal(beta_sd(1, 1), sqrt(1 / 12))
  expect_equal(beta_sd(6.6, 13.4),
               sqrt(6.6 * 13.4 / (20^2 * 21)))
})

test_that("cfr_stan_init only initialises the parameters that exist", {
  base <- list(cfr_a = 2, cfr_b = 3, p1_prior_mean = 2.4, p2_prior_mean = 0.5,
               q1_prior_mean = 3, q2_prior_mean = 0.2)

  # estimated delay, no recovery
  sd1 <- c(base, list(p1_est = 1L, p2_est = 1L, use_recovery = 0L,
                      q1_est = 1L, q2_est = 1L))
  i1 <- cfr_stan_init(sd1)()
  expect_setequal(names(i1), c("cfr", "p1_par", "p2_par"))
  expect_true(i1$cfr > 0 && i1$cfr < 1)
  expect_true(i1$p1_par > 0 && i1$p2_par > 0)

  # fixed delay: only cfr is sampled
  sd2 <- c(base, list(p1_est = 0L, p2_est = 0L, use_recovery = 0L,
                      q1_est = 0L, q2_est = 0L))
  expect_equal(names(cfr_stan_init(sd2)()), "cfr")

  # two-outcome: recovery parameters appear
  sd3 <- c(base, list(p1_est = 1L, p2_est = 1L, use_recovery = 1L,
                      q1_est = 1L, q2_est = 1L))
  expect_setequal(names(cfr_stan_init(sd3)()),
                  c("cfr", "p1_par", "p2_par", "q1_par", "q2_par"))
})

test_that("simulate_linelist requires a delay", {
  expect_error(simulate_linelist(n = 5), "supply a `delay`")
})

test_that("prepare_cfr_data requires onset_date and death_date columns", {
  expect_error(prepare_cfr_data(data.frame(death_date = as.Date("2026-01-05"))),
               "onset_date")
  expect_error(prepare_cfr_data(data.frame(onset_date = as.Date("2026-01-01"))),
               "death_date")
})

test_that("stan_delay_fields rejects non-distspec input and non-native params", {
  expect_error(stan_delay_fields(5, "dist_id", "p"), "distspec distribution")
  # a malformed dist_spec whose parameters are not the family's native ones
  bad <- structure(
    list(parameters = list(foo = 1), distribution = "gamma"),
    class = c("dist_spec", "list")
  )
  expect_error(stan_delay_fields(bad, "dist_id", "p"), "native parameters")
})

test_that("summary.cfrnow_fit rejects objects not from fit_cfr", {
  expect_error(summary.cfrnow_fit(1), "fit_cfr")
})
