#' Perform a chain exchange in parallel tempering
#'
#' This function attempts a swap between two adjacent chains in a parallel tempering
#' setup. It ensures that swaps follow the convention where the higher chain (lower
#' implausibility target) attempts to swap with the lower chain (higher implausibility target).
#'
#' @param X A matrix where each row corresponds to the current position of a chain and each column represents
#' the dimensions of the input space.
#'   The last row represents the chain with the smallest implausibility target.
#' @param imp_levels A list where each element corresponds to the implausibility
#'   levels for each chain (excluding the first). Each element is a vector of
#'   implausibility thresholds for the constraints or waves.
#' @param imp_Xs A matrix of implausibility values where each row corresponds to the implausibility
#' value in the corresponding row of X. Each column represents the implausibility of a wave/constraint. Use
#'   `Inf` for implausibilities that fail earlier constraints.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{Xt}{The updated matrix of chain positions.}
#'   \item{Imp_t}{The updated matrix of implausibility values.}
#' }
#' @details
#' Swaps are always attempted in the downward direction, from a higher chain
#' (higher/easier implausibility target) to a lower chain (lower/nested implausibility target).
#' This is because lower chains automatically satisfy targets on higher chains due to nesting.
#' Randomization is introduced by sampling whether to move `chain_i` up or not.
#' If the lower chain satisfies the implausibility criteria for the higher chain's
#' constraints, the swap occurs.
#'
#' @examples
#' # Example setup
#' # Assume a 2 wave history match with implausibility targets of 3 in each wave. If an input fails
#' # the first wave constraint, its implausibility for the second is `Inf`
#' X <- matrix(rnorm(20), nrow = 5, ncol = 4)
#' print(X)
#' imp_levels <- list(
#'   c(10, Inf),   # Second chain (first has no target), only the first constraint needed
#'   c(3, 8),      # Third chain passes the first constraint so is targeting the space with
#'                 # the first constraint satisfied and the wave 2 implausibility less than 8
#'   c(3, 5),      # Fourth chain almost at the target
#'   c(3, 3)       # The target space
#' )
#'
#' # Generate implausibilities where the second value is `Inf` if the first fails
#' imp_Xs <- rbind(c(9.9, Inf), c(8, Inf), c(2.6, 7.8), c(2.4, 2.5), c(2.2, 2.9))
#' print(imp_Xs)
#'
#' # Perform a chain exchange
#' result <- chain_exchange(X, imp_levels, imp_Xs)
#' print(result)

#' @export
chain_exchange <- function(X, imp_levels, imp_Xs) {
  # Handle edge case where imp_Xs is a vector
  if (is.null(dim(imp_Xs))) {
    dim(imp_Xs) <- c(length(imp_Xs), 1)
  }

  num_chains <- nrow(X)

  # Select a random chain
  chain_i <- sample.int(num_chains, 1)

  if (!(chain_i %in% c(1, num_chains))) {  # Interior chains
    if (runif(1) < 0.5) {  # Random decision
      chain_i <- chain_i + 1  # Move to higher index
    }
    chain_j <- chain_i - 1  # Always the chain directly below
  } else if (chain_i < 2) {  # First chain (boundary case)
    chain_i <- 2
    chain_j <- 1
  } else {  # Last chain (boundary case)
    chain_j <- num_chains - 1
  }

  # Check if the lower chain satisfies the implausibility levels of the higher chain
  num_relevant_imps <- sum(imp_levels[[chain_i - 1]] < Inf)
  if (all(imp_Xs[chain_j, 1:num_relevant_imps] <= imp_levels[[chain_i - 1]][1:num_relevant_imps])) {
    # Swap the positions and implausibilities
    tmp <- X[chain_i, ]
    tmp_imp <- imp_Xs[chain_i, ]
    X[chain_i, ] <- X[chain_j, ]
    imp_Xs[chain_i, ] <- imp_Xs[chain_j, ]
    X[chain_j, ] <- tmp
    imp_Xs[chain_j, ] <- tmp_imp
  }

  list(Xt = X, Imp_t = imp_Xs)
}
