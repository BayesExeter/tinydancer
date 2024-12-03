#' Uniform Slice Sampler over a Box-Region
#'
#' Performs uniform slice sampling over a box-region, avoiding Gibbs steps for computational efficiency in multiple dimensions.
#' This function generates a series of uniform samples using the slice sampling algorithm.
#'
#' @param x0 A numeric vector representing the starting point of the Markov chain.
#' @param f A function to evaluate the target distribution. The function should return -Inf for values outside the feasible region.
#' @param nsmp The number of samples to generate.
#' @param w A numeric vector or scalar of step sizes for each dimension.
#' @param box_limits A matrix with two columns representing the lower and upper bounds for each dimension of the box region.
#' @param ... Additional arguments passed to `f`.
#'
#' @return A matrix where each row is a sample from the target density, with the number of rows equal to `nsmp`.
#' @examples
#' # Example of using slice_sampler_core (USE IMPLAUSIBILITY)
#' f <- function(x) -sum(x^2)  # A simple quadratic function (Gaussian target)
#' x0 <- c(0, 0)
#' nsmp <- 1000
#' w <- c(1, 1)
#' box_limits <- matrix(c(-5, -5, 5, 5), ncol=2)
#' samples <- slice_sampler_core(x0, f, nsmp, w, box_limits)
#' @importFrom stats quantile
#' @importFrom parallel mclapply

#' @export
slice_sampler_core <- function(x0, f, nsmp, w, box_limits, ...) {
  n <- length(x0)
  samples <- matrix(0, nrow=nsmp, ncol=n)

  # Obtain the first sample
  result <- slice_sampler_core_once(x0=x0, f=f, w=w, box_limits=box_limits, ...)
  samples[1,] <- result$sample
  w <- result$w #Use updated w from initial call.

  # Obtain subsequent samples
  for (i in 2:nsmp) {
    result <- slice_sampler_core_once(x0=as.vector(samples[i-1,]), f=f, w=w, box_limits=box_limits, ...)
    samples[i,] <- result$sample
    #w <- result$w #Not used as an initial update of w is normally enough. If a pathological case arises, could explore adapting this and even increasing.
  }

  return(samples)
}

#' Single Step of Uniform Slice Sampling
#'
#' Performs a single step of uniform slice sampling. Given the current position `x0`, the function
#' samples from the target distribution by adaptively adjusting the bounds on each dimension.
#'
#' @param x0 A numeric vector representing the current position in the parameter space.
#' @param f A function to evaluate the target distribution. The function should return -Inf for values outside the feasible region. -Inf is compared to a constant, but any f provided will lead to generation of uniform samples from the region where f > -Inf.
#' @param w A numeric vector or scalar of step sizes for each dimension.
#' @param box_limits A matrix with two columns representing the lower and upper bounds for each dimension of the box region.
#' @param max_retries The maximum number of retries to attempt if the sample fails to meet the constraints. Default is 100.
#' @param ... Additional arguments passed to `f`.
#'
#' @return A list containing:
#'   - `sample`: The sampled point.
#'   - `w`: The adapted step sizes after the sampling step.
#' @note If the sampler exceeds `max_retries`, a warning is issued, and the initial point `x0` is returned as the sample. This ensures continuity but may indicate a problem with the target density function or parameterization.
#' @export
slice_sampler_core_once <- function(x0, f, w, box_limits, max_retries = 100, ...) {
  n <- length(x0)
  if (length(w) == 1)
    w <- rep(w, n)  # Ensure w matches the dimensionality of x0

  # Generate random box around x0 within which to look for new samples
  us <- runif(n)
  ls <- x0 - w * us
  rs <- ls + w

  retry_count <- 0
  while (retry_count < max_retries) {
    #Generate candidate sample
    u_tmp <- runif(n)
    x1 <- ls + u_tmp * (rs - ls)

    # Check if the sampled point is within the box limits
    if (all(x1 >= box_limits[, 1] & x1 <= box_limits[, 2])) {
      # Call f only for feasible x1
      if (f(x1, ...) > -Inf) {
        return(list(sample = x1, w = w))  # Successful sample
      }

      # Update bounds based on density failure
      ls[x1 < x0] <- x1[x1 < x0]
      rs[x1 >= x0] <- x1[x1 >= x0]
    } else {
      # Update bounds based on feasibility failure
      ls[x1 < x0] <- x1[x1 < x0]
      rs[x1 >= x0] <- x1[x1 >= x0]
    }

    # Adaptive w scaling
    if (retry_count %% 10 == 0) {
      w <- 0.9 * w  # Reduce step size after several retries
    }

    retry_count <- retry_count + 1
  }

  warning("Max retries reached and no sample found in slice_sampler_once")
  return(list(sample = x0, w = w))  # Fallback to initial point if max retries reached
}

#' Calculate Log Density
#'
#' This function calculates the log density of a point `x` based on the given
#' implausibility function, levels, and box limits.
#'
#' @param x A numeric vector representing the point for which the log density is to be calculated.
#' @param this_levels A numeric vector specifying the thresholds for implausibility at the current level. `x` will be valid if implausibility is below `this_levels`
#' @param final_levels A numeric vector specifying the target thresholds for implausibility at the final level. Required to pass to implausibility.
#' @param implausibility A user-defined function that calculates the implausibility of a point `x` given the `final_levels`.
#'                       It must accept `x` and `final_levels` as arguments and return a numeric vector.
#'
#' @return Returns `1` if all calculated implausibilities are less than or equal to `this_levels`,
#'         and `-Inf` otherwise.
#'
#' @examples
#' # Define a simple implausibility function
#' implausibility <- function(x, final_levels) { abs(x) - final_levels }
#'
#' # Define parameters
#' x <- c(0.2, 0.3)
#' this_levels <- c(0.3, 0.7)
#' final_levels <- c(0.3, 0.7)
#'
#' # Calculate log density
#' result <- log_dens(x, this_levels, final_levels, implausibility)
#' print(result)
#' @export
log_dens <- function(x, this_levels, final_levels, implausibility) {
  timp <- implausibility(x, final_levels)
  ifelse(all(timp <= this_levels), 1, -Inf)
}

