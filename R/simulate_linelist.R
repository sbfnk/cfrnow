#' Simulate a line list for testing and examples
#'
#' Draws onset dates over a window, marks each case fatal with probability
#' `cfr`, and gives fatal cases an onset-to-death delay drawn from `delay`. Pass
#' a `recovery` delay to also simulate onset-to-recovery times for the non-fatal
#' cases and add a `recovery_date` column. Delays are dist.spec distributions
#' with fixed parameters, matching what [fit_cfr()] takes; you can specify them
#' by mean and sd (e.g. `dist.spec::LogNormal(mean = 12.75, sd = 7)`). Returns a
#' line list with `onset_date`, `death_date` (`NA` for non-fatal cases) and,
#' when `recovery` is given, `recovery_date` (`NA` for fatal cases). The full,
#' untruncated outcomes are simulated; pass the result to [prepare_cfr_data()]
#' with an `obs_time` to induce the real-time truncation.
#'
#' @param n Number of cases.
#' @param cfr True case fatality ratio.
#' @param delay Onset-to-death delay: a dist.spec distribution
#'   ([dist.spec::LogNormal()] or [dist.spec::Gamma()]) with fixed parameters.
#' @param recovery Optional onset-to-recovery delay (same form as `delay`); when
#'   given, non-fatal cases get a `recovery_date`.
#' @param onset_start First possible onset date.
#' @param onset_days Width of the onset window (days); onsets are uniform
#'   over it.
#' @return A data frame with `onset_date`, `death_date` and, if `recovery` is
#'   given, `recovery_date`.
#' @examples
#' simulate_linelist(n = 5, cfr = 0.6,
#'                   delay = dist.spec::LogNormal(mean = 12.75, sd = 7))
#' @export
simulate_linelist <- function(n = 200, cfr = 0.5, delay, recovery = NULL,
                              onset_start = as.Date("2026-01-01"),
                              onset_days = 60) {
  if (missing(delay)) {
    stop("supply a `delay` (a dist.spec distribution with fixed parameters).",
         call. = FALSE)
  }
  onset <- as.Date(onset_start) + sample.int(onset_days, n, replace = TRUE) - 1
  fatal <- stats::runif(n) < cfr

  otd <- sample_delay(n, delay)
  death_date <- as.Date(rep(NA, n))
  death_date[fatal] <- onset[fatal] + round(otd[fatal])
  out <- data.frame(onset_date = onset, death_date = death_date)

  if (!is.null(recovery)) {
    otr <- sample_delay(n, recovery)
    recovery_date <- as.Date(rep(NA, n))
    recovery_date[!fatal] <- onset[!fatal] + round(otr[!fatal])
    out$recovery_date <- recovery_date
  }
  out
}
