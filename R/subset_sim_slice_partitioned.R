#' Perform subset simulation with slice sampling on partitioned space.
#'
#' This function divides the sampling space into partitions, performs subset
#' simulation with slice sampling within each partition, and combines results
#' from successful partitions.
#'
#' @param implausibility A function to calculate implausibility for a given point.
#' @param dims Integer. Dimensionality of the sampling space.
#' @param target_levels Numeric. Target implausibility level(s) to achieve.
#' @param control_list List. Control parameters for the subset simulation. Default options include:
#'   \itemize{
#'     \item `num_switches`: Number of switches in the input space.
#'     \item `volume_ratio`: Volume ratio for sampling.
#'     \item `num_mutations`: Number of mutations per iteration.
#'     \item `max_num_chains`: Maximum number of chains to use per partition.
#'     \item `box_limits`: Matrix defining the bounds of the sampling space.
#'     \item `switch_settings_list`: List of settings for switches.
#'     \item `debug_mode`: Logical. If TRUE, run in debug mode.
#'     \item `levels_dp`: Decimal places to round implausibility levels.
#'     \item `one_per_level`: Logical. If TRUE, only retain one sample per level.
#'   }
#' @param n_partitions Integer. Number of partitions to divide the sampling space into.
#' @return A list containing combined samples and their corresponding implausibilities.
#' @examples
#' # Example usage
#' st1 <- rbind(c(0.000002, 0.000000875), c(0.000000875, 0.00025))
#' st2 <- rbind(c(0.00005, 0.000000825), c(0.000000825, 0.000002))
#' sdtiny <- function(x, m1 = c(0, 0), m2 = c(5, 4),
#'                    s1 = st1, s2 = st2) {
#'   coef1 <- (1 / ((s1[2, 2] * s1[1, 1]) - (s1[2, 1] * s1[1, 2])))
#'   s1inv <-  coef1 * rbind(c(s1[2, 2], -s1[1, 2]), c(-s1[2, 1], s1[1, 1]))
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
#' result <- subset_sim_slice_partitioned(
#'   implausibility,
#'   dims = 2,
#'   target_levels = 3,
#'   control_list = control_list,
#'   n_partitions = 2
#' )
#' @details
#' Use this function to generate uniform samples within partitions. To combine into uniform samples of the target space, an importance sampling wrapper is required that estimates the subvolumes of each partition containing samples. *Due for next release*
#'
#' @export
subset_sim_slice_partitioned <- function(implausibility, dims, target_levels = 3, control_list = list(), n_partitions = 4) {
  # Define defaults
  defaults <- list(
    num_switches = 0,
    volume_ratio = 0.1,
    num_mutations = 30,
    max_num_chains = 8,
    box_limits = NULL,
    switch_settings_list = NULL,
    debug_mode = FALSE,
    levels_dp = 2,
    one_per_level = FALSE
  )

  # Merge defaults with user-specified options
  control_list <- modifyList(defaults, control_list)
  box_limits <- control_list$box_limits

  # Partition space
  partitions <- partition_space(box_limits, n_partitions)

  # Run simulations on partitions
  results <- if (control_list$debug_mode) {
    lapply(partitions, function(partition) {
      control_list$box_limits <- partition
      subset_sim_slice(implausibility, dims, target_levels, control_list)
    })
  } else {
    mclapply(partitions, function(partition) {
      control_list$box_limits <- partition
      subset_sim_slice(implausibility, dims, target_levels, control_list)
    })
  }

  # Filter successful results
  successful_results <- Filter(function(res) res$reached_target, results)
  if (length(successful_results) == 0) {
    stop("No partitions successfully reached the target! Consider allowing more partitions or more chains per partition.")
  }

  # Combine results from all partitions
  combined_samples <- do.call(rbind, lapply(successful_results, function(res) res$X))
  #combined_implausibilities <- do.call(c, lapply(successful_results, function(res) res$Implausibilities))
  combined_implausibilities <- do.call(
    if (is.vector(successful_results[[1]]$Implausibilities)) c else rbind,
    lapply(successful_results, function(res) res$Implausibilities)
  )

  return(list(X = combined_samples, Implausibilities = combined_implausibilities))
}


#' Partition the Sampling Space
#'
#' This function partitions a multidimensional box into smaller subspaces. Each partition is represented
#' as a set of lower and upper bounds for each dimension.
#'
#' @param box_limits Numeric matrix. A `dims x 2` matrix where each row specifies the lower and upper bounds
#'        of a dimension in the sampling space.
#' @param n_partitions Integer. The number of partitions to create along each dimension.
#'
#' @return A list where each element is a matrix representing a partition. Each matrix has two columns,
#'         `"Lower"` and `"Upper"`, and rows corresponding to dimensions.
partition_space <- function(box_limits, n_partitions) {
  dims <- nrow(box_limits)

  # Generate partition edges for each dimension
  partition_edges <- lapply(1:dims, function(d) {
    seq(box_limits[d, 1], box_limits[d, 2], length.out = n_partitions + 1)
  })

  # Generate all combinations of indices for partitions
  partition_indices <- expand.grid(rep(list(1:n_partitions), dims))

  # Construct partitions from indices
  partitions <- lapply(1:nrow(partition_indices), function(i) {
    indices <- as.numeric(partition_indices[i, ])
    bounds <- matrix(0, nrow = dims, ncol = 2)
    for (d in 1:dims) {
      bounds[d, 1] <- partition_edges[[d]][indices[d]]       # Lower bound
      bounds[d, 2] <- partition_edges[[d]][indices[d] + 1]   # Upper bound
    }
    colnames(bounds) <- c("Lower", "Upper")
    return(bounds)
  })

  return(partitions)
}
