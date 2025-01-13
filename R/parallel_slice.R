#' Parallel Tempering Slice Sampling for Uniform target regions
#'
#' This function performs parallel tempering MCMC with slice sampling for uniform samples from a level set of a function (called implausibility throughout).
#'
#' @param num_chains Integer. The number of parallel chains.
#' @param num_mutations Integer. The number of slice samples taken per chain per iteration (the term mutation arises from the Evolutionary Monte Carlo literature).
#' @param num_iterations Integer. The total number of iterations.
#' @param d Integer. The dimensionality of the parameter space.
#' @param imp_levels List of numeric vectors. The implausibility levels for each chain (excluding the first chain), see `note` below.
#' @param final_target_levels Numeric vector. The final implausibility levels for all constraints.
#' @param x_starts Matrix. `num_chains` x `d`. Starting positions for all chains except the first.
#' @param box_limits Matrix. The box region where sampling is constrained; each row corresponds to the lower and upper bounds of a dimension.
#' @param debug_mode Logical. If `TRUE`, runs in debug mode with `lapply`; otherwise uses `mclapply` for parallelism. MAY NEED TO CHANGE FOR WINDOWS
#' @param num_switches Integer. Number of switches for discrete settings in the first chain (default: `0`).
#' @param switch_settings_list List. A list of POSSIBLE settings for switches in the first chain.
#' @param volume_ratio Numeric. Ratio for scaling the step size `w` across chains.
#' @param print_every Integer. Frequency of status updates during iterations (default: `100`).
#' @param implausibility A user-defined function that calculates the implausibility
#' of a point given the required levels. (see main examples)
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{uniform_sample_list}}{A list of matrices, one per chain, containing the sampled values. The last chain contains samples from the `final_target_levels`. The matrices are `num_iterations`x `d+num_waves` where `num_waves` is the length of an element of `imp_levels`.}
#'   \item{\code{restart}}{The state of the chains for restarting the sampler.}
#' }
#' @note The first chain samples uniformly over the box, so it does not use implausibility levels or starting values.
#'
parallel_slice <- function(num_chains, num_mutations, num_iterations, d, imp_levels,
                           final_target_levels, x_starts, box_limits, debug_mode = FALSE,
                           num_switches = 0, switch_settings_list = NULL, volume_ratio, print_every = 100,
                           implausibility) {
  # Determine the number of waves
  num_waves <- length(imp_levels[[length(imp_levels)]])

  # Preallocate chains
  chain_list <- vector("list", num_chains)
  for (zz in 1:num_chains) chain_list[[zz]] <- matrix(0, nrow = num_iterations, ncol = d + num_waves)

  if (!is.matrix(x_starts)) dim(x_starts) <- c(1, d)

  # Initial uniform sample for the first chain
  chain1_sample <- sapply(1:nrow(box_limits), function(i) runif(1, min = box_limits[i, 1], max = box_limits[i, 2]))
  if (num_switches > 0) {
    if (is.null(switch_settings_list))
      stop("Please provide a switch settings list.")
    switches <- unlist(lapply(switch_settings_list, function(e) sample(e, 1)))
    chain1_sample <- c(chain1_sample, switches)
  }
  chain1_imp <- implausibility(chain1_sample, final_target_levels)

  # Parallel sampling for other chains
  parallel_chains <- if (debug_mode) {
    lapply(2:num_chains, function(j) one_slice(
      x_start = x_starts[j - 1, ],
      M = num_mutations,
      w = volume_ratio^((j - 1) / d),
      box_limits = box_limits,
      imp_level = imp_levels[[j - 1]],
      final_levels = final_target_levels,
      final_num_waves = num_waves,
      implausibility = implausibility
    ))
  } else {
    plan("multicore", wokers=num_chains-1)
    future.apply::future_lapply(2:num_chains, function(j) one_slice(
      x_start = x_starts[j - 1, ],
      M = num_mutations,
      w = volume_ratio^((j - 1) / d),
      box_limits = box_limits,
      imp_level = imp_levels[[j - 1]],
      final_levels = final_target_levels,
      final_num_waves = num_waves,
      implausibility = implausibility
    ), future.seed = TRUE)
  }

  # Preallocate combined chain results
  xn_imp <- matrix(0, nrow = num_chains, ncol = d + num_waves)
  xn_imp[1, ] <- c(chain1_sample, chain1_imp)
  for (k in 2:num_chains) {
    xn_imp[k, ] <- c(parallel_chains[[k - 1]]$Xt, parallel_chains[[k - 1]]$Imp_t)
  }

  chain_info <- list(Xt = xn_imp[, 1:d], Imp_t = xn_imp[, (d + 1):(d + num_waves)])

  # Exchange step
  num_swaps <- max(1, num_chains - 1)
  for (i in seq_len(num_swaps)) {
    chain_info <- chain_exchange(chain_info$Xt, imp_levels, chain_info$Imp_t)
  }

  # Store the initial results
  for (a in seq_len(num_chains)) {
    chain_list[[a]][1, ] <- c(chain_info$Xt[a, ], chain_info$Imp_t[a, ])
  }

  # Main iteration loop
  for (k in seq(2, num_iterations)) {
    # Mutation step for first chain
    chain1_sample <- sapply(1:nrow(box_limits), function(i) runif(1, min = box_limits[i, 1], max = box_limits[i, 2]))
    if (num_switches > 0) {
      switches <- unlist(lapply(switch_settings_list, function(e) sample(e, 1)))
      chain1_sample <- c(chain1_sample, switches)
    }
    chain1_imp <- implausibility(chain1_sample, final_target_levels)

    # Parallel mutation step
    parallel_chains <- if (debug_mode) {
      lapply(2:num_chains, function(j) one_slice(
        x_start = as.vector(chain_info$Xt[j, ]),
        M = num_mutations,
        w = volume_ratio^((j - 1) / d),
        box_limits = box_limits,
        imp_level = imp_levels[[j - 1]],
        final_levels = final_target_levels,
        final_num_waves = num_waves,
        implausibility = implausibility
      ))
    } else {
      plan("multicore", workers = num_chains-1)
      future.apply::future_lapply(2:num_chains, function(j) one_slice(
        x_start = as.vector(chain_info$Xt[j, ]),
        M = num_mutations,
        w = volume_ratio^((j - 1) / d),
        box_limits = box_limits,
        imp_level = imp_levels[[j - 1]],
        final_levels = final_target_levels,
        final_num_waves = num_waves,
        implausibility = implausibility
      ), future.seed = TRUE)
    }

    xn_imp[1, ] <- c(chain1_sample, chain1_imp)
    for (l in seq(2, num_chains)) {
      xn_imp[l, ] <- c(parallel_chains[[l - 1]]$Xt, parallel_chains[[l - 1]]$Imp_t)
    }

    chain_info <- list(Xt = xn_imp[, 1:d], Imp_t = xn_imp[, (d + 1):(d + num_waves)])

    # Exchange step
    for (i in seq_len(num_swaps)) {
      chain_info <- chain_exchange(chain_info$Xt, imp_levels, chain_info$Imp_t)
    }

    # Store the results
    for (a in seq_len(num_chains)) {
      chain_list[[a]][k, ] <- c(chain_info$Xt[a, ], chain_info$Imp_t[a, ])
    }

    if (k %% print_every == 0) {
      message(sprintf("Iteration %d of %d complete.", k, num_iterations))
    }
  }

  return(list(uniform_sample_list = chain_list, restart = parallel_chains))
}



#' One Slice Sampling Step for a Single Chain for parallel tempering
#'
#' This function performs slice sampling for a single chain within a parallel tempering algorithm using the `slice_sampler_core` function. (Often called mutation steps in the evolutionary Monte Carlo literature).
#'
#' @param M Integer. The number of samples to draw.
#' @param imp_level Numeric vector. The implausibility levels for this chain.
#' @param final_levels Numeric vector. The final implausibility levels for all constraints.
#' @param x_start Numeric vector. The starting position for the chain.
#' @param box_limits Matrix. The box region where sampling is constrained; each row corresponds to a dimension.
#' @param final_num_waves Integer. The total number of implausibility constraints (waves).
#' @param w Numeric or numeric vector. The initial step size(s) for the slice sampler.
#' @param implausibility A user-defined function that calculates the implausibility
#' of a point given the required levels. (see main examples)
#' @param ... Additional arguments to pass to the log-density function `f`.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{Xt}}{Numeric vector. The final sampled position after the mutation step.}
#'   \item{\code{Imp_t}}{Numeric vector. The implausibility values for the sampled position.}
#' }
#' @note This function assumes that `f`, the log-density function, is defined globally or passed via \code{...}.
#' It also assumes that the implausibility function \code{implausibility} is globally defined and evaluates the constraints at the given position.
one_slice <- function(M, imp_level, final_levels, x_start, box_limits, final_num_waves, w = 1, implausibility, ...) {
  # Perform slice sampling for M steps
  sampled_points <- slice_sampler_core(
    x0 = x_start,
    f = log_dens, # Assumes log_dens is defined globally or passed via ...
    nsmp = M,
    w = w,
    box_limits = box_limits,
    this_levels = imp_level,
    final_levels = final_levels,
    implausibility = implausibility,
    ...
  )

  # Extract the final sample and compute its implausibility
  final_sample <- sampled_points[M, ]
  final_imp <- implausibility(final_sample, final_levels) # Assumes implausibility is globally defined

  # Return results
  return(list(
    Xt = final_sample,
    Imp_t = final_imp
  ))
}
