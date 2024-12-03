#' Perform subset simulation with slice sampling
#'
#' This function initializes and performs Markov Chain Monte Carlo (MCMC) sampling using a
#' subset simulation approach to identify implausibility levels and generate samples.
#'
#' @param implausibility A function that takes a vector of length \code{dims} as the first argument
#' and computes the function for which there is a level set target. Output is a single value or vector of values the same length as the length of `target_levels`
#' @param dims Integer. The dimensionality of the parameter space.
#' @param target_levels Numeric vector. A vector of implausibility levels describing the cutoff for each wave.
#' A vector of length 1 indicates a single wave and could be used to indicate a single level set of a function.
#' @param control_list List. A list of control parameters with the following elements:
#'        \itemize{
#'          \item `num_switches`: Integer. Number of switch inputs (default: 0).
#'          \item `volume_ratio`: Numeric. Volume ratio for slice sampling (default: 0.1).
#'          \item `num_mutations`: Integer. Number of mutations per chain (default: 30).
#'          \item `num_iterations`: Integer. Number of iterations for the final sampling wave (default: 1000).
#'          \item `box_limits`: Numeric matrix (dims x 2). Bounds for the initial space.
#'          \item `switch_settings_list`: List of length `num_switches` containing the possible values for each switch input (default: NULL).
#'          \item `debug_mode`: Logical. Enable debug messages (default: FALSE).
#'          \item `max_num_chains`: Integer. Maximum number of chains (default: 8).
#'          \item `levels_dp`: Integer. Decimal places for implausibility levels (default: 2).
#'          \item `one_per_level`: Logical. Whether to generate a single chain targeting each value of `target_levels` separately (default: FALSE).
#'        }
#' @return A list containing:
#'         \itemize{
#'           \item `x`: Numeric matrix. Uniform samples from the final sampling wave.
#'           \item `implausibilities`: Numeric vector or matrix. Implausibility values associated with the samples.
#'           \item `reached_target`: Logical. Whether the target implausibility levels were reached.
#'         }
#' @examples
#' # Example usage
#' st1 = rbind(c(0.000002, 0.000000875), c(0.000000875, 0.00025))
#' st2 = rbind(c(0.00005, 0.000000825), c(0.000000825, 0.000002))
#' sdtiny <- function(x, m1 = c(0, 0), m2 = c(5, 4),
#'                    s1 = st1, s2 = st2) {
#'   coef1 <- (1 / ((s1[2, 2] * s1[1, 1]) - (s1[2, 1] * s1[1, 2])))
#'   s1inv <- coef1 * rbind(c(s1[2, 2], -s1[1, 2]), c(-s1[2, 1], s1[1, 1]))
#'   coef2 <- (1 / ((s2[2, 2] * s2[1, 1]) - (s2[2, 1] * s2[1, 2])))
#'   s2inv <- coef2 * rbind(c(s2[2, 2], -s2[1, 2]), c(-s2[2, 1], s2[1, 1]))
#'   sdevs1 <- sqrt(t(x - m1) %*% s1inv %*% (x - m1))
#'   sdevs2 <- sqrt(t(x - m2) %*% s2inv %*% (x - m2))
#'   min(sdevs1, sdevs2)
#' }
#' implausibility <- function(x, target_level = 3) {
#'   sdtiny(x)
#' }
#' control_list <- list(
#'   box_limits = cbind(rep(-3, 2), rep(7, 2)),
#'   num_mutations = 40,
#'   num_iterations = 100
#' )
#' result <- subset_sim_slice(
#'   implausibility,
#'   dims = 2,
#'   target_levels = 3,
#'   control_list = control_list
#' )

#' @details
#' This function is very effective if confident that the target space is simply connected. If not, and one or more disconnected regions is suspected, the `subset_sim_slice_partitioned` function will enable additional regions to be found if using a number of partitions. If used on a simply connected region, this function will return uniform samples.
#'
#' @export
subset_sim_slice <- function(implausibility, dims, target_levels = 3, control_list = list()) {
  # Set default control parameters
  defaults <- list(
    num_switches = 0,
    volume_ratio = 0.1,
    num_mutations = 30,
    num_iterations = 1000,
    box_limits = NULL,
    switch_settings_list = NULL,
    debug_mode = FALSE,
    max_num_chains = 8,
    levels_dp = 2,
    one_per_level = FALSE
  )

  control_list <- modifyList(defaults, control_list)

  # Check box_limits
  if (!dim(control_list$box_limits)[1] == dims - control_list$num_switches) {
    stop("box_limits must be a dims x 2 matrix describing the support of the initial space")
  } else {
    box_limits <- control_list$box_limits
  }

  # Initialize first chain
  chain1_sample <- sapply(1:(dims - control_list$num_switches), function(i)
    runif(control_list$num_mutations * 10, min = box_limits[i, 1], max = box_limits[i, 2]))
  num_waves <- length(target_levels)
  num_switches <- control_list$num_switches
  switch_settings_list <- control_list$switch_settings_list

  if (num_switches > 0) {
    switches <- t(sapply(1:(control_list$num_mutations * 10), function(i)
      unlist(lapply(switch_settings_list, function(e) sample(e, 1)))))
    chain1_sample <- cbind(chain1_sample, switches)
  }

  # Compute implausibilities for the first chain
  chain1_imps <- t(sapply(1:nrow(chain1_sample), function(k)
    implausibility(chain1_sample[k, ], target_levels)))

  if (num_waves == 1) {
    x_starts <- chain1_sample[which.min(chain1_imps), 1:dims]
    dim(x_starts) <- c(1, dims)
    current_lower <- max(c(target_levels, min(chain1_imps)))
  } else {
    # Multi-wave initialization
    t_indices <- 1:control_list$num_mutations
    in_space_i <- TRUE
    i <- 0
    while (in_space_i & (i < num_waves)) {
      i <- i + 1
      in_space_i <- any(chain1_imps[t_indices, i] <= target_levels[i])
    }
    if (i == 1) {
      start_index <- which.min(chain1_imps[, 1])
      current_lower <- chain1_imps[start_index, ]
    } else if (i < num_waves) {
      current_lower <- target_levels
      start_index <- which.min(chain1_imps[, i])
      current_lower[i:num_waves] <- chain1_imps[start_index, i:num_waves]
    } else {
      current_lower <- target_levels
      start_index <- which.min(chain1_imps[, i])
      if (!in_space_i) {
        current_lower[i] <- chain1_imps[start_index, i]
      }
    }
    x_starts <- chain1_sample[start_index, ]
    dim(x_starts) <- c(1, dims)
  }

  # Main loop
  num_chains <- 1
  reached_target <- all(current_lower <= target_levels)
  while (any(current_lower > target_levels) & num_chains <= control_list$max_num_chains) {
    if (control_list$debug_mode) {
      print(paste("Current number of chains = ", num_chains, ". Seeking the next implausibility level...", sep = ""))
    }

    new_samples <- one_slice_all(
      m = control_list$num_mutations,
      imp_level = current_lower,
      final_levels = target_levels,
      x_start = x_starts[, 1:dims],
      box_limits = box_limits,
      final_num_waves = num_waves,
      w = control_list$volume_ratio^(num_chains / dims),
      implausibility = implausibility
    )

    # Update levels and starting points
    if (num_waves == 1) {
      current_lower <- max(c(target_levels, min(new_samples[, dims + 1])))
      x_starts <- new_samples[which.min(new_samples[, dims + 1]), 1:dims]
      dim(x_starts) <- c(1, dims)
    } else {
      t_indices <- 1:control_list$num_mutations
      in_space_i <- TRUE
      i <- 0
      while (in_space_i & (i < num_waves)) {
        i <- i + 1
        in_space_i <- any(new_samples[t_indices, dims + i] <= target_levels[i])
      }
      if (i == 1) {
        start_index <- which.min(new_samples[, dims + 1])
        current_lower <- new_samples[start_index, -c(1:dims)]
      } else if (i < num_waves) {
        current_lower <- target_levels
        start_index <- which.min(new_samples[, dims + i])
        current_lower[i:num_waves] <- new_samples[start_index, (dims + i):(dims + num_waves)]
      } else {
        current_lower <- target_levels
        start_index <- which.min(new_samples[, dims + i])
        if (!in_space_i) {
          current_lower[i] <- new_samples[start_index, dims + i]
        }
      }
      x_starts <- new_samples[start_index, 1:dims]
      dim(x_starts) <- c(1, dims)
    }

    if (control_list$debug_mode) {
      print(paste("This level implausibility found = ", list(current_lower), sep = " "))
    }
    num_chains <- num_chains + 1
    reached_target <- all(current_lower <= target_levels)
  }

  if (reached_target) {
    if (control_list$debug_mode) {
      print("Generating Final Samples")
    }
    final_samples <- one_slice_all(
      m = control_list$num_iterations,
      imp_level = current_lower,
      final_levels = target_levels,
      x_start = x_starts[, 1:dims],
      box_limits = box_limits,
      final_num_waves = num_waves,
      w = control_list$volume_ratio^(num_chains / dims),
      implausibility = implausibility
    )
    return(list(
      X = final_samples[, 1:dims],
      Implausibilities = final_samples[, -c(1:dims)],
      reached_target = TRUE
    ))
  } else {
    return(list(
      X = x_starts[, 1:dims],
      Implausibilities = current_lower,
      reached_target = FALSE
    ))
  }
}
