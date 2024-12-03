#' Add Rung to Temperature Ladder for Optimal Mixing in Parallel Tempering MCMC
#'
#' This function adds a temperature/implausibility level to a current Parallel Tempering MCMC that has not yet warmed up to the target level.
#'
#' @param chain Matrix. The chain containing samples at the previous lowest level of a PTMCMC step (see `parallel_slice`). Contains `d` + `length(target_imps)` columns, where the first `d` entries in a row are the sample and the final entries are the target function (e.g. Implausibility) values at the different waves of the corresponding sample.
#' @param target_imps Numeric vector. Target implausibility levels for each chain.
#' @param d Integer. Dimensionality of the parameter space.
#' @param p Numeric. The quantile probability used to set the next level set (implausibility level).
#' @param levels_dp Integer. The number of decimal places for rounding the level set evaluation.
#' @param num_at_target Integer. The number of chains already fixed at the target level set values.
#' @param one_per_level Logical. If `TRUE`, ensures a unique chain for each target implausibility level. Set to `FALSE` for a more efficient search for the final target level in all waves and to `TRUE` if the goal is to have uniform samples for each wave for plotting or analysis. Defaults to `FALSE`.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{new_imps}}{The adjusted implausibility levels.}
#'   \item{\code{start_value}}{The starting value for the next iteration.}
#'   \item{\code{num_at_target}}{The updated number of chains fixed at the target implausibility levels.}
#' }
#'
#' @examples
#' # Example usage:
#' set.seed(123)
#' chain <- matrix(runif(100), nrow = 10, ncol = 10)
#' target_imps <- c(0.1, 0.5, 1, 2, 5, Inf)
#' d <- 3
#' p <- 0.5
#' levels_dp <- 2
#' num_at_target <- 2
#' one_per_level <- TRUE
#' result <- add_temperature_rung(
#'   chain = chain,
#'   target_imps = target_imps,
#'   d = d,
#'   p = p,
#'   levels_dp = levels_dp,
#'   num_at_target = num_at_target,
#'   one_per_level = one_per_level
#' )
#' print(result)
#' @export
add_temperature_rung <- function(chain, target_imps, d, p, levels_dp,
                                       num_at_target, one_per_level = FALSE) {
  # Determine new implausibility levels for chains beyond those fixed at target
  new_imps <- sapply((num_at_target + 1):length(target_imps), function(i) {
    round(quantile(chain[, d + i], probs = c(0, p, 1), names = FALSE)[2], digits = levels_dp)
  })

  new_num_at_target <- num_at_target

  # Ensure new implausibility levels do not fall below target levels
  for (i in (num_at_target + 1):length(target_imps)) {
    if (new_imps[i - num_at_target] < target_imps[i]) {
      new_imps[i - num_at_target] <- target_imps[i]
      new_num_at_target <- new_num_at_target + 1
      # Enforce uniqueness of chains per target implausibility level
      if (one_per_level) break
    }
  }

  # Add fixed target implausibility levels if any
  if (num_at_target > 0) {
    new_imps <- c(target_imps[1:num_at_target], new_imps)
  }

  # Set remaining implausibility levels to infinity if applicable
  if (new_num_at_target < (length(target_imps) - 1)) {
    new_imps[(new_num_at_target + 2):length(target_imps)] <- Inf
  }

  # Determine starting value based on the sampled value with the lowest implausibility level
  if (new_num_at_target < 1) {
    start_value <- chain[which.min(chain[, d + 1]), ]
  } else if (new_num_at_target == length(target_imps)) {
    start_value <- chain[which.min(chain[, d + length(target_imps)]), ]
  } else {
    start_value <- chain[which.min(chain[, d + 1 + new_num_at_target]), ]
  }

  return(list(
    new_imps = new_imps,
    start_value = start_value,
    num_at_target = new_num_at_target
  ))
}
