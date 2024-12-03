#' One Slice All: Sampling from a Slice Sampler and Returning All Samples
#'
#' This function performs slice sampling and returns all the samples along with their corresponding implausibility levels.
#' It uses `slice_sampler_core` for sampling and `Implausibility` to calculate implausibility values for each sample.
#'
#' @param m Integer. The number of samples to generate.
#' @param imp_level Numeric vector. The current implausibility level for sampling.
#' @param final_levels Numeric vector. The final target implausibility levels.
#' @param x_start Numeric vector. The starting point for the slice sampler.
#' @param box_limits Numeric matrix. The limits of the sampling box.
#' @param final_num_waves Integer. The number of final waves to consider in the sampling process.
#' @param w Numeric. The width parameter used in the slice sampling process. Default is 1.
#' @param implausibility A user-defined function that calculates the implausibility
#' of a point given the required levels. (see main examples)
#' @return A matrix where each row corresponds to a sample, and columns represent the sampled values along with
#'         the corresponding implausibility levels.
#' @examples
#' # Example implausibility function
#' implausibility <- function(x, final_levels) {
#'   abs(x) # Placeholder: Euclidean norm as implausibility
#' }
#' result <- one_slice_all(
#'   m = 1000,
#'   imp_level = 0.5,
#'   final_levels = c(0.3, 0.7),
#'   x_start = c(0.2, 0.3),
#'   box_limits = matrix(c(0, 1, 0, 1), ncol = 2),
#'   final_num_waves = 2,
#'   implausibility=implausibility
#' )
#' @usage one_slice_all(m, imp_level, final_levels, x_start,
#'  box_limits, final_num_waves, w = 1, implausibility)
#' @export
one_slice_all <- function(m, imp_level, final_levels, x_start, box_limits, final_num_waves, w = 1, implausibility) {
  # Perform the slice sampling
  t_slice <- slice_sampler_core(x0 = x_start, f = log_dens, nsmp = m, w = w, box_limits = box_limits,
                                this_levels = imp_level, final_levels = final_levels, implausibility=implausibility)

  # Calculate implausibility for each sample
  timps <- sapply(1:m, function(k) implausibility(t_slice[k, ], final_levels))

  # Return samples with implausibility values
  if (length(final_levels) > 1) {
    chain_sample <- cbind(t_slice, t(timps))
  } else {
    chain_sample <- cbind(t_slice, timps)
  }

  return(chain_sample)
}
