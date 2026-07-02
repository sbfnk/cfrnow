#' Supported onset-to-death delay families
#'
#' The families `cfrnow` can fit, each an analytical primarycensored delay
#' distribution paired with a uniform primary event. These match the dist.spec
#' constructors accepted by [fit_cfr()]'s `delay` argument
#' ([dist.spec::LogNormal()], [dist.spec::Gamma()]).
#'
#' @return A character vector of supported family names.
#' @examples
#' cfrnow_families()
#' @export
cfrnow_families <- function() {
  c("lognormal", "gamma")
}

# Match a family name against the supported set (single source of validation).
validate_family <- function(delay_family) {
  match.arg(delay_family, cfrnow_families())
}

# Draw n delays (days) from a family parameterised by its mean and sd, used by
# simulate_linelist() for both onset-to-death and onset-to-recovery delays.
sample_delay <- function(n, mean, sd, delay_family) {
  native <- delay_to_native(mean, sd, delay_family)
  switch(
    validate_family(delay_family),
    lognormal = stats::rlnorm(n, native[["meanlog"]], native[["sdlog"]]),
    gamma = stats::rgamma(n, shape = native[["shape"]], rate = native[["rate"]])
  )
}

# Convert a delay (mean, sd) to a family's native parameters, for
# simulate_linelist() so simulations match the family the model fits.
delay_to_native <- function(mean, sd, delay_family) {
  delay_family <- validate_family(delay_family)
  switch(
    delay_family,
    lognormal = {
      cv2 <- (sd / mean)^2
      c(meanlog = log(mean) - 0.5 * log1p(cv2), sdlog = sqrt(log1p(cv2)))
    },
    gamma = c(shape = (mean / sd)^2, rate = mean / sd^2)
  )
}
