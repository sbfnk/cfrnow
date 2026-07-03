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

# Draw n delays (days) from a dist.spec distribution with fixed parameters, used
# by simulate_linelist() for the onset-to-death and onset-to-recovery delays.
sample_delay <- function(n, delay) {
  fam <- dist.spec::get_distribution(delay)
  pars <- dist.spec::get_parameters(delay)[delay_native_order(fam)]
  if (!all(vapply(pars, is.numeric, logical(1)))) {
    stop("simulate_linelist() needs a delay with fixed parameters (numbers), ",
         "not priors.", call. = FALSE)
  }
  switch(
    fam,
    lognormal = stats::rlnorm(n, pars[["meanlog"]], pars[["sdlog"]]),
    gamma = stats::rgamma(n, shape = pars[["shape"]], rate = pars[["rate"]])
  )
}
