# Vendor the primarycensored Stan functions into inst/stan/include/.
#
# Run this whenever primarycensored is updated. It writes a self-contained
# functions file (with dependencies resolved) that inst/stan/curecfr.stan
# #includes, so the package does not depend on the user's installed
# primarycensored library layout at compile time. Commit the generated file.

# functions = NULL vendors the whole primarycensored function set. We only call
# primarycensored_lpmf() and primarycensored_cdf(), but the latter references
# primarycensored_ode by name in its numerical-integration branch (unused for
# analytical lognormal/gamma/weibull + uniform, but it must still be in scope at
# compile time). Vendoring everything keeps the include self-contained.
primarycensored::pcd_load_stan_functions(
  functions = NULL,
  wrap_in_block = FALSE,
  write_to_file = TRUE,
  output_file = file.path("inst", "stan", "include", "pcd_functions.stan")
)
