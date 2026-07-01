#' Simulate a line list for testing and examples
#'
#' Draws onset dates over a window, marks each case fatal with probability
#' `cfr`, and gives fatal cases an onset-to-death delay from the chosen family
#' (parameterised by its mean and sd, matching the model). Returns a line list
#' with `onset_date` and `death_date` (`NA` for non-fatal cases). The full,
#' untruncated outcomes are simulated; pass the result to [prepare_cfr_data()]
#' with an `obs_time` to induce the real-time truncation.
#'
#' @param n Number of cases.
#' @param cfr True case fatality ratio.
#' @param delay_mean,delay_sd Onset-to-death delay mean and sd (days).
#' @param delay_family Delay family; see [curecfr_families()].
#' @param onset_start First possible onset date.
#' @param onset_days Width of the onset window (days); onsets are uniform over it.
#' @return A data frame with `onset_date` and `death_date`.
#' @examples
#' simulate_linelist(n = 5, cfr = 0.6)
#' @export
simulate_linelist <- function(n = 200, cfr = 0.5,
                              delay_mean = 12.75, delay_sd = 7,
                              delay_family = "gamma",
                              onset_start = as.Date("2026-01-01"),
                              onset_days = 60) {
  native <- delay_to_native(delay_mean, delay_sd, delay_family)
  onset <- as.Date(onset_start) + sample.int(onset_days, n, replace = TRUE) - 1
  fatal <- stats::runif(n) < cfr
  otd <- switch(
    match.arg(delay_family, curecfr_families()),
    lognormal = stats::rlnorm(n, native[["meanlog"]], native[["sdlog"]]),
    gamma = stats::rgamma(n, shape = native[["shape"]], rate = native[["rate"]])
  )
  death_date <- as.Date(rep(NA_real_, n), origin = "1970-01-01")
  death_date[fatal] <- onset[fatal] + round(otd[fatal])
  data.frame(onset_date = onset, death_date = death_date)
}
