#' Native parameter order for a delay family
#'
#' The two native parameters in the positional order used by the delay
#' distribution. Doubles as cfrnow's supported-family guard (lognormal, gamma).
#' @param family Delay family name.
#' @return Character vector of the two native parameter names.
#' @noRd
delay_native_order <- function(family) {
  native <- switch(family,
    lognormal = c("meanlog", "sdlog"),
    gamma = c("shape", "rate")
  )
  if (is.null(native)) {
    stop("unsupported delay family '", family,
      "'; cfrnow supports lognormal and gamma.",
      call. = FALSE
    )
  }
  native
}
