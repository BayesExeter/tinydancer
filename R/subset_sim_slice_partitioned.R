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
#'     \item `num_iterations`: Integer. Number of iterations for the final sampling wave (default: 1000).
#'     \item `max_num_chains`: Maximum number of chains to use per partition.
#'     \item `box_limits`: Matrix defining the bounds of the sampling space.
#'     \item `switch_settings_list`: List of settings for switches.
#'     \item `debug_mode`: Logical. If TRUE, run each subset_sim in debug mode (with progress messages).
#'     \item `forking`: Logical. If TRUE, each subset simulation is run in parallel with mclapply. Defaults to FALSE, see `details` section below.
#'     \item `levels_dp`: Decimal places to round implausibility levels.
#'     \item `one_per_level`: Logical. If TRUE, only retain one sample per level.
#'   }
#' @param n_partitions Integer. Number of partitions to divide the sampling space into.
#' @param random Logical. If TRUE, dimensions are chosen randomly for splitting. Default is FALSE.
#' @return A list containing combined samples, `X`, and their corresponding implausibilities, `Implausibilities`. The list also includes an estimate of the relative volume of the target space, `volume_estimate`.
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
#' Uniform samples are generated within each partition using slice sampling. Note uniformity is guaranteed up to the guarantees afforded by slice samplers (we may not have converged/found all of the target compatible space). Samples are then combined across partitions using importance sampling with weights estimated according to the relative volumes of the subspaces found within each partition, to return uniform samples from the target space.
#'
#' If 2 or more partitions are requested, slice sampling will try to run in parallel via mclapply if `control_list$fork_chains` is `TRUE` and via `lapply` otherwise (the default). Forked processing in R is not stable, not available to all operating systems (e.g. it does not work on windows) and, at the time of writing, Rstudio does not permit it via the future package. Depending on the complexity of `implausibility`, forking can be highly effective, particularly when the function is simple. However, functions with calls to Rcpp or Python backends may not work and `mclapply` may not return any errors or may even hang.
#'
#' @export
subset_sim_slice_partitioned <- function(implausibility, dims, target_levels = 3, control_list = list(), n_partitions = 4, random=FALSE) {
  # Define defaults
  defaults <- list(
    num_switches = 0,
    volume_ratio = 0.1,
    num_mutations = 30,
    num_iterations = 1000,
    max_num_chains = 8,
    box_limits = NULL,
    estimate_volume = TRUE,
    switch_settings_list = NULL,
    debug_mode = FALSE,
    forking = FALSE,
    levels_dp = 2,
    one_per_level = FALSE
  )

  # Merge defaults with user-specified options
  control_list <- modifyList(defaults, control_list)
  box_limits <- control_list$box_limits
  if(!control_list$estimate_volume){
    control_list$estimate_volume <- TRUE
    message("estimate_volume is required for uniform sampling and has been turned on")
  }

  # Partition space
  partitions <- partition_space(box_limits, n_partitions)
  partition_volumes <- lapply(partitions, function(e) prod(e[,'Upper'] - e[,'Lower']))

  # Run simulations on partitions
  results <- if (!control_list$forking) {
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

  for(jj in seq_along(results)){
    results[[jj]]$target_volume <- results[[jj]]$volume_estimate*partition_volumes[[jj]]
  }

  # Filter successful results
  successful_results <- Filter(function(res) res$reached_target, results)
  if (length(successful_results) == 0) {
    stop("No partitions successfully reached the target! Consider allowing more partitions or more chains per partition.")
  }

  #Importance sampling to return a uniform sample
  target_volume <- sum(sapply(successful_results, function(kk) kk$target_volume))
  importance_weights <- sapply(successful_results, function(res) res$target_volume)/target_volume
  sampled_list_indices <- sample(seq_along(successful_results),
                                 size = control_list$num_iterations, replace = TRUE, prob= importance_weights)
  counts <- table(sampled_list_indices)
  names(counts) <- as.numeric(names(counts))

  # Sampling from the lists
  sampled_results <- lapply(names(counts), function(idx) {
    idx <- as.numeric(idx)  # Convert name back to numeric index
    list_element <- successful_results[[idx]]  # Get the corresponding list element

    # Number of samples to take
    n_samples <- counts[[as.character(idx)]]

    # Sample from rows of X and corresponding Implausibilities
    n_rows <- if (is.matrix(list_element$X)) nrow(list_element$X) else length(list_element$X)
    sampled_rows <- sample(seq_len(n_rows), size = n_samples, replace = FALSE)

    sampled_X <- list_element$X[sampled_rows, , drop = FALSE]

    if (is.vector(list_element$Implausibilities)) {
      # If Implausibilities is a vector, extract sampled elements directly
      sampled_Implausibilities <- list_element$Implausibilities[sampled_rows]
    } else {
      # If Implausibilities is a matrix, extract sampled rows
      sampled_Implausibilities <- list_element$Implausibilities[sampled_rows, , drop = FALSE]
    }

    list(X = sampled_X, Implausibilities = sampled_Implausibilities)
  })

  # Combine results into a single structure
  combined_samples <- do.call(rbind, lapply(sampled_results, function(res) res$X))
  combined_implausibilities <- do.call(
    if (is.vector(successful_results[[1]]$Implausibilities)) c else rbind,
    lapply(sampled_results, function(res) res$Implausibilities)
  )

  #Produce a relative volume (of the target to the full space) estimate
  relative_volume <- target_volume/prod(box_limits[,2] - box_limits[,1])

  return(list(X = combined_samples, Implausibilities = combined_implausibilities, volume_estimate = relative_volume))
}


  # Combine results from all partitions
#  combined_samples <- do.call(rbind, lapply(successful_results, function(res) res$X))
#  combined_implausibilities <- do.call(
#    if (is.vector(successful_results[[1]]$Implausibilities)) c else rbind,
#    lapply(successful_results, function(res) res$Implausibilities)
#  )

#  return(list(X = combined_samples, Implausibilities = combined_implausibilities))
#}


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
