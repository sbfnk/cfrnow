#' Simulate a line list for testing and examples
#'
#' Draws onset dates over a window, marks each case fatal with probability
#' `cfr`, and gives fatal cases an onset-to-death delay from the chosen family
#' (parameterised by its mean and sd, matching the model). With
#' `recovery = TRUE` the non-fatal cases also get an onset-to-recovery delay and
#' a `recovery_date` column; recoveries resolve more slowly than deaths by
#' default,
#' which is what drives the resolved-set enrichment the model corrects for.
#' Returns a line list with `onset_date`, `death_date` (`NA` for non-fatal
#' cases) and, when requested, `recovery_date` (`NA` for fatal cases). The full,
#' untruncated outcomes are simulated; pass the result to [prepare_cfr_data()]
#' with an `obs_time` to induce the real-time truncation.
#'
#' @param n Number of cases.
#' @param cfr True case fatality ratio.
#' @param delay_mean,delay_sd Onset-to-death delay mean and sd (days).
#' @param delay_family Onset-to-death delay family; see [cfrnow_families()].
#' @param recovery If `TRUE`, also simulate non-fatal recoveries and add a
#'   `recovery_date` column.
#' @param recovery_mean,recovery_sd Onset-to-recovery delay mean and sd (days),
#'   used when `recovery = TRUE`.
#' @param recovery_family Onset-to-recovery delay family; see
#'   [cfrnow_families()].
#' @param onset_start First possible onset date.
#' @param onset_days Width of the onset window (days); onsets are uniform
#'   over it.
#' @return A data frame with `onset_date`, `death_date` and, if
#'   `recovery = TRUE`, `recovery_date`.
#' @examples
#' simulate_linelist(n = 5, cfr = 0.6)
#' simulate_linelist(n = 5, cfr = 0.6, recovery = TRUE)
#' @export
simulate_linelist <- function(n = 200, cfr = 0.5,
                              delay_mean = 12.75, delay_sd = 7,
                              delay_family = "gamma",
                              recovery = FALSE,
                              recovery_mean = 21, recovery_sd = 9,
                              recovery_family = "gamma",
                              onset_start = as.Date("2026-01-01"),
                              onset_days = 60) {
  onset <- as.Date(onset_start) + sample.int(onset_days, n, replace = TRUE) - 1
  fatal <- stats::runif(n) < cfr

  otd <- sample_delay(n, delay_mean, delay_sd, delay_family)
  death_date <- as.Date(rep(NA, n))
  death_date[fatal] <- onset[fatal] + round(otd[fatal])
  out <- data.frame(onset_date = onset, death_date = death_date)

  if (recovery) {
    otr <- sample_delay(n, recovery_mean, recovery_sd, recovery_family)
    recovery_date <- as.Date(rep(NA, n))
    recovery_date[!fatal] <- onset[!fatal] + round(otr[!fatal])
    out$recovery_date <- recovery_date
  }
  out
}
