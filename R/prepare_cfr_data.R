#' Build model inputs from a line list
#'
#' Classifies each case at an observation cut-off into an observed death (with
#' an interval-censored onset-to-death delay), a right-censored survivor, or a
#' fully-resolved non-death, and returns the pieces [fit_cfr()] passes to Stan.
#'
#' Onset is taken over the day-window `[onset_lower, onset_upper]` when those
#' columns are present (defaulting to a one-day window at `onset_date`). Deaths
#' are recorded to the day.
#'
#' @param linelist A data frame with an `onset_date` column, an optional
#'   `onset_lower`/`onset_upper` onset window, and a `death_date` column
#'   (`NA` for cases that have not died). Dates may be `Date` or coercible.
#' @param obs_time Real-time cut-off (`Date` or coercible), or `NULL` for a
#'   retrospective fit in which every recorded death counts and survivors are
#'   treated as fully resolved. In real time, deaths dated after `obs_time` are
#'   treated as not-yet-known and the case is right-censored.
#' @param t0 Optional time origin (`Date`). Defaults to `min(onset) - 60` days.
#' @param max_delay Largest plausible onset-to-death delay (days). Death records
#'   implying a longer or negative delay are dropped as date-entry errors.
#'
#' @return A `curecfr_data` list with the Stan inputs (`n_death`,
#'   `death_delay`, `death_width`, `n_cens`, `censor_time`, `censor_width`,
#'   `n_resolved`) plus `n_cases`, `n_deaths`, `t0` and `obs_time`.
#' @examples
#' ll <- simulate_linelist(n = 50)
#' prepare_cfr_data(ll, obs_time = as.Date("2026-02-01"))
#' @export
prepare_cfr_data <- function(linelist, obs_time = NULL, t0 = NULL,
                             max_delay = 60) {
  if (!"onset_date" %in% names(linelist)) {
    stop("`linelist` needs an `onset_date` column.", call. = FALSE)
  }
  if (!"death_date" %in% names(linelist)) {
    stop("`linelist` needs a `death_date` column (NA for non-fatal cases).",
         call. = FALSE)
  }

  onset <- as.Date(linelist$onset_date)
  onset_lo <- if ("onset_lower" %in% names(linelist)) {
    as.Date(linelist$onset_lower)
  } else {
    onset
  }
  onset_up <- if ("onset_upper" %in% names(linelist)) {
    as.Date(linelist$onset_upper)
  } else {
    onset
  }
  death <- as.Date(linelist$death_date)

  if (is.null(t0)) t0 <- min(onset, na.rm = TRUE) - 60
  t0 <- as.Date(t0)

  onset_lo_day <- as.numeric(onset_lo - t0)
  width <- as.numeric(onset_up - onset_lo) + 1        # >= 1 day
  retrospective <- is.null(obs_time)
  obs_offset <- if (retrospective) Inf else as.numeric(as.Date(obs_time) - t0)

  death_day <- as.numeric(death - t0)                 # NA for non-fatal
  is_death <- !is.na(death_day) & death_day <= obs_offset
  delay <- death_day - onset_lo_day                   # from onset-window start

  # Drop impossible onset->death delays (date-entry errors) among the deaths.
  bad <- is_death & (delay < 0 | delay > max_delay)
  n_dropped <- sum(bad, na.rm = TRUE)
  if (n_dropped > 0) {
    warning(sprintf("dropped %d onset->death date-entry error(s)", n_dropped),
            call. = FALSE)
  }
  keep <- !bad
  is_death <- is_death & keep
  is_surv <- keep & !is_death

  death_delay <- as.integer(round(delay[is_death]))
  death_width <- as.numeric(width[is_death])

  if (retrospective) {
    n_resolved <- sum(is_surv)
    censor_time <- numeric(0)
    censor_width <- numeric(0)
  } else {
    n_resolved <- 0L
    censor_time <- pmax(obs_offset - onset_lo_day[is_surv], 0)
    censor_width <- as.numeric(width[is_surv])
  }

  structure(
    list(
      n_death = length(death_delay),
      death_delay = death_delay,
      death_width = death_width,
      n_cens = length(censor_time),
      censor_time = as.numeric(censor_time),
      censor_width = censor_width,
      n_resolved = as.integer(n_resolved),
      n_cases = sum(keep),
      n_deaths = length(death_delay),
      t0 = t0,
      obs_time = if (retrospective) NA else as.Date(obs_time)
    ),
    class = "curecfr_data"
  )
}
