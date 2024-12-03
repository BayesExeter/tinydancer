#' Perform Parallel Tempering MCMC Slice Sampling for uniform target regions
#'
#' This function performs uniform sampling in the specified space using Parallel Tempering MCMC with slice sampling.
#' It utilizes the temperature ladder constructed by the provided sampler object.
#'
#' @param sampler A list returned by `construct_temperature_ladder`, containing the implausibility levels
#'        (`imp_levels`), starting points (`x_starts`), and control parameters (`control_list`).
#' @param n_iter Integer. Number of iterations for the sampler. Default is 1000.
#' @param implausibility A user-defined function that calculates the implausibility
#' of a point given the required levels. (see main examples)
#' @param control_list List. A list of control parameters that can override the settings in `sampler$control_list`.
#'        Available options include:
#'        \itemize{
#'          \item `num_mutations`: Number of mutations per iteration.
#'          \item `debug_mode`: Logical, whether to enable debug mode.
#'          \item `print_every`: Number of iterations between progress messages.
#'        }
#'
#' @return A list containing samples from the Parallel Tempering MCMC slice sampler.
#' @examples
#' # Example usage
#' st1 <- rbind(c(0.000002, 0.000000875), c(0.000000875, 0.00025))
#' st2 <- rbind(c(0.00005, 0.000000825), c(0.000000825, 0.000002))
#' sdtiny <- function(x, m1=c(0,0), m2=c(5,4),
#'                    s1=st1, s2=st2) {
#'   coef1 <- (1 / ((s1[2, 2] * s1[1, 1]) - (s1[2, 1] * s1[1, 2])))
#'   s1inv <- coef1 * rbind(c(s1[2, 2], -s1[1, 2]), c(-s1[2, 1], s1[1, 1]))
#'   coef2 <- (1 / ((s2[2, 2] * s2[1, 1]) - (s2[2, 1] * s2[1, 2])))
#'   s2inv <- coef2 * rbind(c(s2[2, 2], -s2[1, 2]), c(-s2[2, 1], s2[1, 1]))
#'   sdevs1 <- sqrt(t(x - m1) %*% s1inv %*% (x - m1))
#'   sdevs2 <- sqrt(t(x - m2) %*% s2inv %*% (x - m2))
#'   return(min(sdevs1, sdevs2))
#' }
#' implausibility <- function(x, target_level=3) {
#'   sdtiny(x)
#' }
#' control_list <- list(
#'   box_limits = cbind(rep(-3, 2), rep(7, 2)),
#'   num_mutations = 8,
#'   num_iterations = 100
#' )
#' new_ladder <- construct_temperature_ladder(
#'   implausibility = implausibility,
#'   dims = 2,
#'   target_levels = 3,
#'   control_list = control_list
#' )
#' result <- sample_slice_ptmcmc(new_ladder, n_iter = 1000, implausibility=implausibility)
#' @export
sample_slice_ptmcmc <- function(sampler, n_iter = 1000, control_list = list(), implausibility) {
  # Merge default control values with user-specified options
  control_list$num_mutations <- control_list$num_mutations %||% sampler$control_list$num_mutations
  control_list$debug_mode <- control_list$debug_mode %||% sampler$control_list$debug_mode
  control_list$print_every <- control_list$print_every %||% sampler$control_list$print_every

  # Retrieve dimensions and temperature ladder details
  n_waves <- length(sampler$imp_levels[[1]])
  d <- ncol(sampler$x_starts) - n_waves

  # Perform sampling using the Parallel Slice method
  parallel_slice(
    num_chains = length(sampler$imp_levels) + 1,
    num_mutations = control_list$num_mutations,
    num_iterations = n_iter,
    d = d,
    imp_levels = sampler$imp_levels,
    final_target_levels = sampler$control_list$target_level,
    x_starts = sampler$x_starts[, -c((d + 1):(d + n_waves))],
    box_limits = sampler$control_list$box_limits,
    debug_mode = control_list$debug_mode,
    num_switches = sampler$control_list$num_switches,
    switch_settings_list = sampler$control_list$switch_settings_list,
    volume_ratio = sampler$control_list$volume_ratio,
    print_every = control_list$print_every,
    implausibility=implausibility
  )
}


`%||%` <- function(a, b) {
  if (!is.null(a)) a else b
}
