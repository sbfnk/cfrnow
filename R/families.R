#' Supported onset-to-death delay families
#'
#' The families `cfrnow` can fit, each an analytical primarycensored delay
#' distribution paired with a uniform primary event and parameterised here by
#' its mean and standard deviation. Weibull is analytically available in
#' primarycensored but not yet wired into the mean/sd conversion, so it is
#' excluded here.
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

# Validate a family name and return its primarycensored delay dist_id.
delay_dist_id <- function(delay_family) {
  primarycensored::pcd_stan_dist_id(
    validate_family(delay_family), type = "delay"
  )
}

# Convert a delay (mean, sd) to a family's native (p1, p2), mirroring the Stan
# `delay_to_native()`. Used by simulate_linelist() so simulations match the
# family the model fits.
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
