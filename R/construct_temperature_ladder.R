#' Construct Temperature Ladder for Parallel Tempering MCMC
#'
#' This function initializes the full Markov Chain (parallel chains) using
#' slice sampling and parallel tempering. The temperature ladder is constructed iteratively
#' by finding implausibility levels for each wave that ensure good mixing.
#'
#' @param implausibility A function that takes a vector of length \code{dims} as the first argument
#' and computes the function for which there is a level set target. Output is a single value or vector of values the same length as the length of `target_levels`.
#' @param dims Integer. Dimensionality of the parameter space.
#' @param target_levels Numeric vector. A vector of implausibility levels describing the cutoff for each wave.
#' A vector of length 1 indicates a single wave and could be used to indicate a single level set of a function.
#' @param control_list List. A list of control parameters with defaults:
#' \describe{
#'   \item{\code{num_switches}}{Integer. Number of switch variables in the model. Default: 0.}
#'   \item{\code{volume_ratio}}{Numeric. Ratio for volume-based thresholding. Increase to improve mixing between chains and decrease to use fewer chains to find the target space. Default: 0.1.}
#'   \item{\code{num_mutations}}{Integer. Number of iterations (often called mutations in the Evolutionary Monte Carlo literature) for individual chains within the parallel structure. Default: 30.}
#'   \item{\code{num_iterations}}{Integer. Number of iterations of the paralel chains. Default: 1000.}
#'   \item{\code{box_limits}}{Matrix. A \code{dims x 2} matrix specifying the lower and upper bounds of the initial space.}
#'   \item{\code{switch_settings_list}}{List. A list of length num_switches, each element containing a vector of the possible settings of each switch. Defaults to NULL}
#'   \item{\code{levels_dp}}{Integer. Number of decimal places for implausibility levels. Default: 2.}
#'   \item{\code{one_per_level}}{Logical. If \code{TRUE}, ensures unique implausibility levels for each wave. Default: \code{FALSE}.}
#'   \item{\code{print_every}}{Integer. Frequency of progress updates during sampling to generate new temperature levels. Default: 100.}
#' }
#'
#' @return A list with:
#'   - `imp_levels`: List of temperature/implausibility levels for each wave.
#'   - `x_starts`: Matrix of starting values for the chains.
#'   - `control_list`: Updated control list with additional elements used during sampling.
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
#'
#' # Example with multiple levels of implausibility (e.g., from multiple waves of history matching)
#' ss1 <- rbind(c(0.0002, 0.0000875), c(0.0000875, 0.0025))
#' ss2 <- rbind(c(0.005, 0.0000825), c(0.0000825, 0.0002))
#' sdsmall <- function(x, m1=c(0,0), m2=c(5,4),
#'                     s1=ss1, s2=ss2) {
#'   coef1 <- (1 / ((s1[2, 2] * s1[1, 1]) - (s1[2, 1] * s1[1, 2])))
#'   s1inv <- coef1 * rbind(c(s1[2, 2], -s1[1, 2]), c(-s1[2, 1], s1[1, 1]))
#'   coef2 <- (1 / ((s2[2, 2] * s2[1, 1]) - (s2[2, 1] * s2[1, 2])))
#'   s2inv <- coef2 * rbind(c(s2[2, 2], -s2[1, 2]), c(-s2[2, 1], s2[1, 1]))
#'   sdevs1 <- sqrt(t(x - m1) %*% s1inv %*% (x - m1))
#'   sdevs2 <- sqrt(t(x - m2) %*% s2inv %*% (x - m2))
#'   return(min(sdevs1, sdevs2))
#' }
#' sb1 <- rbind(c(0.002, 0.000875), c(0.000875, 0.025))
#' sb2 <- rbind(c(0.05, 0.000825), c(0.000825, 0.002))
#' sdbig <- function(x, m1=c(0,0), m2=c(5,4),
#'                   s1=sb1, s2=sb2) {
#'   coef1 <- (1 / ((s1[2, 2] * s1[1, 1]) - (s1[2, 1] * s1[1, 2])))
#'   s1inv <- coef1 * rbind(c(s1[2, 2], -s1[1, 2]), c(-s1[2, 1], s1[1, 1]))
#'   coef2 <- (1 / ((s2[2, 2] * s2[1, 1]) - (s2[2, 1] * s2[1, 2])))
#'   s2inv <- coef2 * rbind(c(s2[2, 2], -s2[1, 2]), c(-s2[2, 1], s2[1, 1]))
#'   sdevs1 <- sqrt(t(x - m1) %*% s1inv %*% (x - m1))
#'   sdevs2 <- sqrt(t(x - m2) %*% s2inv %*% (x - m2))
#'   return(min(sdevs1, sdevs2))
#' }
#' implausibility <- function(x, targetLevel, levels=5, waves=3) {
#'   ans <- rep(Inf, levels)
#'   waveFail <- FALSE
#'   this.level <- 1
#'   wave.num <- 1
#'   Timp <- NA
#'   while ((this.level <= levels) & !waveFail) {
#'     if (wave.num == 1) {
#'       Timp <- sdbig(x)
#'     } else if (wave.num == 2) {
#'       Timp <- sdsmall(x)
#'     } else {
#'       Timp <- sdtiny(x)
#'     }
#'     wave.num <- wave.num + 1
#'     if ((Timp > targetLevel[this.level]) & (wave.num <= waves)) {
#'       waveFail <- TRUE
#'     }
#'     if ((!waveFail) & (wave.num > waves)) {
#'       ans[this.level:levels] <- Timp
#'       this.level <- levels + 1
#'     } else {
#'       ans[this.level] <- Timp
#'       this.level <- this.level + 1
#'     }
#'   }
#'   return(ans)
#' }
#' control_list <- list(
#'   num_mutations = 8,
#'   num_iterations = 100,
#'   box_limits = cbind(rep(-3, 2), rep(7, 2))
#' )
#' new_ladder <- construct_temperature_ladder(
#'   implausibility = implausibility,
#'   dims = 2,
#'   target_levels = c(3, 3, 3, 2.5, 2),
#'   control_list = control_list
#' )
#'@details Currently no stopping rule exists for this, though the printing is sufficiently verbose for the user to understand what is happening. If a target compatible subspace exists, this algorithm will find it eventually, however, if the ladder seems to be converging far above the target, it may be worth the user intervening.
#'@export
construct_temperature_ladder <- function(implausibility, dims, target_levels, control_list = list()) {
  # Define defaults
  defaults <- list(
    num_switches = 0,
    volume_ratio = 0.1,
    num_mutations = 30,
    num_iterations = 1000,
    box_limits = NULL,
    switch_settings_list = NULL,
    debug_mode = FALSE,
    levels_dp = 2,
    one_per_level = FALSE,
    print_every = 100
  )

  # Merge user-supplied `control_list` with defaults
  control_list <- modifyList(defaults, control_list)

  # Validate box limits
  if (!is.matrix(control_list$box_limits) || nrow(control_list$box_limits) != (dims - control_list$num_switches)) {
    stop("`box_limits` must be a (dims - num_switches) x 2 matrix describing the support of the initial space.")
  }
  box_limits <- control_list$box_limits

  num_waves <- length(target_levels) # Determine the number of waves
  num_iterations <- control_list$num_iterations
  num_switches <- control_list$num_switches
  switch_settings_list <- control_list$switch_settings_list

  # Initialize chain 1
  chain1_sample <- sapply(1:(dims - num_switches), function(i)
    runif(num_iterations, min = box_limits[i, 1], max = box_limits[i, 2]))

  if (num_switches > 0) {
    if (is.null(switch_settings_list)) {
      stop("`switch_settings_list` must be provided if `num_switches` > 0.")
    }
    switches <- t(sapply(1:num_iterations, function(i)
      unlist(lapply(switch_settings_list, function(e) sample(e, 1)))))
    chain1_sample <- cbind(chain1_sample, switches)
  }

  # Compute implausibilities
  chain1_imps <- sapply(1:num_iterations, function(k)
    implausibility(chain1_sample[k, ], target_levels))

  if (num_waves > 1) {
    chain1_sample <- cbind(chain1_sample, t(chain1_imps))
  } else {
    chain1_sample <- cbind(chain1_sample, chain1_imps)
  }

  # Initialize lists and variables
  all_chains <- list(chain1_sample)
  new_level <- add_temperature_rung(chain = chain1_sample, target_imps = target_levels,
                                    d = dims, p = control_list$volume_ratio,
                                    levels_dp = control_list$levels_dp,
                                    num_at_target = 0,
                                    one_per_level = control_list$one_per_level)
  current_lower <- new_level$new_imps
  imp_levels <- list(current_lower)
  if(length(target_levels)<2)
    message(sprintf("Currently sampling a ladder with level sets: %.2f. Seeking the next ladder rung...", imp_levels))
  else
    message(sprintf("Currently sampling a ladder with level sets: %s. Seeking the next ladder rung...", paste(sprintf("%.2f", unlist(imp_levels)), collapse = ", ")))
  x_starts <- new_level$start_value
  dim(x_starts) <- c(1, dims + num_waves)
  num_chains <- 2
  debug_mode <- control_list$debug_mode

  # Construct the temperature ladder
  while (any(current_lower > target_levels)) {
    message(sprintf("Current number of chains = %d. Seeking the next implausibility level...", num_chains))
    run_chains <- parallel_slice(
      num_chains = num_chains,
      num_mutations = control_list$num_mutations,
      num_iterations = num_iterations,
      d = dims,
      imp_levels = imp_levels,
      final_target_levels = target_levels,
      x_starts = x_starts[, -c((dims + 1):(dims + num_waves))],
      box_limits = box_limits,
      debug_mode = debug_mode,
      num_switches = num_switches,
      switch_settings_list = switch_settings_list,
      volume_ratio = control_list$volume_ratio,
      print_every = control_list$print_every,
      implausibility = implausibility
    )

    # Update chains and levels
    for (i in 1:(num_chains - 1)) {
      all_chains[[i]] <- rbind(all_chains[[i]], run_chains$uniform_sample_list[[i]])
    }
    all_chains[[num_chains]] <- run_chains$uniform_sample_list[[num_chains]]

    new_level <- add_temperature_rung(chain = all_chains[[num_chains]], target_imps = target_levels,
                                      d = dims, p = control_list$volume_ratio,
                                      levels_dp = control_list$levels_dp,
                                      num_at_target = new_level$num_at_target,
                                      one_per_level = control_list$one_per_level)
    current_lower <- new_level$new_imps
    imp_levels <- c(imp_levels, list(current_lower))
    if(length(target_levels)<2)
      message(sprintf("Currently sampling a ladder with level sets: %s. Seeking the next ladder rung...",
                    paste(sprintf("%.2f", imp_levels), collapse = ", ")))
    else
      message(sprintf(
        "Currently sampling a ladder with level sets: %s. Seeking the next ladder rung...",
        paste(
          lapply(imp_levels, function(levels) paste(sprintf("%.2f", levels), collapse = ", ")),
          collapse = " / "
        )
      ))
      #message(sprintf("Currently sampling a ladder with level sets: %s. Seeking the next ladder rung...", paste(sprintf("%.2f", unlist(imp_levels)), collapse = ", ")))

    x_starts <- matrix(0, nrow = num_chains - 1, ncol = dims + num_waves)
    for (i in 1:(num_chains - 1)) {
      x_starts[i, ] <- all_chains[[i + 1]][nrow(all_chains[[i + 1]]), ]
    }
    x_starts <- rbind(x_starts, new_level$start_value)
    num_chains <- num_chains + 1
  }

  message("Ladder construction complete.")
  message(sprintf("Final temperature ladder: %s", paste(imp_levels, collapse = ", ")))

  return(list(
    imp_levels = imp_levels,
    x_starts = x_starts,
    control_list = c(control_list, list(implausibility = implausibility, target_levels = target_levels))
  ))
}

