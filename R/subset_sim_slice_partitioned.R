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
#' @param random Logical. If TRUE, dimensions are chosen randomly for splitting. Default is FALSE.
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
subset_sim_slice_partitioned <- function(implausibility, dims, target_levels = 3, control_list = list(), n_partitions = 4, random=FALSE) {
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
#' This function partitions a multidimensional box into exactly `n_partitions` subspaces.
#' The splitting prioritizes dimensions sequentially: the first dimension is divided first,
#' followed by the second, and so on, as needed. Alternatively, the splits can be done in a random order.
#'
#' @param box_limits Numeric matrix. A dims x 2 matrix where each row specifies the lower and upper bounds
#'        of a dimension in the sampling space.
#' @param n_partitions Integer. The exact number of partitions to create.
#' @param random Logical. If TRUE, dimensions are chosen randomly for splitting. Default is FALSE.
#'
#' @return A list where each element is a matrix representing a partition. Each matrix has two columns,
#'         "Lower" and "Upper", and rows corresponding to dimensions.
partition_space <- function(box_limits, n_partitions, random = FALSE) {
  dims <- nrow(box_limits)

  # Initialize with the entire box as the first partition
  partitions <- list(box_limits)

  # Track the dimension to split
  split_dim <- 1
  partition_index <- 1

  # Initialize dimension order for random splitting
  if (random) {
    dim_order <- sample(1:dims)
    dim_index <- 1
  }

  # Keep splitting partitions until we reach the desired number
  while (length(partitions) < n_partitions) {
    # Find the partition to split
    current_partition <- partitions[[partition_index]]

    # Choose the dimension to split
    if (random) {
      split_dim <- dim_order[dim_index]
      dim_index <- dim_index + 1
      if (dim_index > dims) {
        dim_order <- sample(1:dims)
        dim_index <- 1
      }
    }

    # Compute the midpoint for the split in the current dimension
    mid_point <- mean(current_partition[split_dim, ])

    # Create two new partitions by splitting along the selected dimension
    partition1 <- current_partition
    partition2 <- current_partition
    partition1[split_dim, 2] <- mid_point  # Update upper bound of the first partition
    partition2[split_dim, 1] <- mid_point  # Update lower bound of the second partition

    # Replace the current partition with the two new partitions
    partitions[[partition_index]] <- partition1
    partitions <- append(partitions, list(partition2))

    # Move to the next dimension for splitting if not random
    if (!random) {
      split_dim <- (split_dim %% dims) + 1
    }

    partition_index <- partition_index + 1
  }

  # Ensure column names are set for all partitions
  partitions <- lapply(partitions, function(part) {
    colnames(part) <- c("Lower", "Upper")
    part
  })

  return(partitions)
}
