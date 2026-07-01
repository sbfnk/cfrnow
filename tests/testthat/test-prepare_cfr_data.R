test_that("retrospective fit resolves every non-death and counts every death", {
  ll <- data.frame(
    onset_date = as.Date("2026-01-01") + c(0, 1, 2, 3),
    death_date = as.Date(c("2026-01-10", NA, "2026-01-15", NA))
  )
  d <- prepare_cfr_data(ll, obs_time = NULL)
  expect_equal(d$n_deaths, 2L)
  expect_equal(d$n_resolved, 2L)      # both survivors treated as resolved
  expect_equal(d$n_cens, 0L)          # nothing right-censored retrospectively
  expect_equal(d$n_cases, 4L)
})

test_that("real-time cut-off censors survivors and hides later deaths", {
  ll <- data.frame(
    onset_date = as.Date("2026-01-01") + c(0, 1, 2),
    # death 1 known by cut-off; death 2 happens after the cut-off; case 3 alive
    death_date = as.Date(c("2026-01-05", "2026-01-20", NA))
  )
  d <- prepare_cfr_data(ll, obs_time = as.Date("2026-01-10"))
  expect_equal(d$n_deaths, 1L)        # only the death dated on/before cut-off
  expect_equal(d$n_cens, 2L)          # the future death + the still-alive case
  expect_equal(d$n_resolved, 0L)
  expect_true(all(d$censor_time >= 0))
})

test_that("impossible onset->death delays are dropped with a warning", {
  ll <- data.frame(
    onset_date = as.Date(c("2026-01-10", "2026-01-01")),
    death_date = as.Date(c("2026-01-05", "2026-01-08"))  # first: death before onset
  )
  expect_warning(d <- prepare_cfr_data(ll, obs_time = NULL), "date-entry")
  expect_equal(d$n_deaths, 1L)
  expect_equal(d$n_cases, 1L)         # the bad record is dropped entirely
})

test_that("onset windows widen the primary censoring width", {
  ll <- data.frame(
    onset_date = as.Date("2026-01-01"),
    onset_lower = as.Date("2026-01-01"),
    onset_upper = as.Date("2026-01-04"),
    death_date = as.Date("2026-01-12")
  )
  d <- prepare_cfr_data(ll, obs_time = NULL)
  expect_equal(d$death_width, 4)      # (upper - lower) + 1
})
