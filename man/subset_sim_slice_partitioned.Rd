% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/subset_sim_slice_partitioned.R
\name{subset_sim_slice_partitioned}
\alias{subset_sim_slice_partitioned}
\title{Perform subset simulation with slice sampling on partitioned space.}
\usage{
subset_sim_slice_partitioned(
  implausibility,
  dims,
  target_levels = 3,
  control_list = list(),
  n_partitions = 4
)
}
\arguments{
\item{implausibility}{A function to calculate implausibility for a given point.}

\item{dims}{Integer. Dimensionality of the sampling space.}

\item{target_levels}{Numeric. Target implausibility level(s) to achieve.}

\item{control_list}{List. Control parameters for the subset simulation. Default options include:
\itemize{
\item \code{num_switches}: Number of switches in the input space.
\item \code{volume_ratio}: Volume ratio for sampling.
\item \code{num_mutations}: Number of mutations per iteration.
\item \code{max_num_chains}: Maximum number of chains to use per partition.
\item \code{box_limits}: Matrix defining the bounds of the sampling space.
\item \code{switch_settings_list}: List of settings for switches.
\item \code{debug_mode}: Logical. If TRUE, run in debug mode.
\item \code{levels_dp}: Decimal places to round implausibility levels.
\item \code{one_per_level}: Logical. If TRUE, only retain one sample per level.
}}

\item{n_partitions}{Integer. Number of partitions to divide the sampling space into.}
}
\value{
A list containing combined samples and their corresponding implausibilities.
}
\description{
This function divides the sampling space into partitions, performs subset
simulation with slice sampling within each partition, and combines results
from successful partitions.
}
\details{
Use this function to generate uniform samples within partitions. To combine into uniform samples of the target space, an importance sampling wrapper is required that estimates the subvolumes of each partition containing samples. \emph{Due for next release}
}
\examples{
# Example usage
st1 <- rbind(c(0.000002, 0.000000875), c(0.000000875, 0.00025))
st2 <- rbind(c(0.00005, 0.000000825), c(0.000000825, 0.000002))
sdtiny <- function(x, m1 = c(0, 0), m2 = c(5, 4),
                   s1 = st1, s2 = st2) {
  coef1 <- (1 / ((s1[2, 2] * s1[1, 1]) - (s1[2, 1] * s1[1, 2])))
  s1inv <-  coef1 * rbind(c(s1[2, 2], -s1[1, 2]), c(-s1[2, 1], s1[1, 1]))
  coef2 <- (1 / ((s2[2, 2] * s2[1, 1]) - (s2[2, 1] * s2[1, 2])))
  s2inv <- coef2 * rbind(c(s2[2, 2], -s2[1, 2]), c(-s2[2, 1], s2[1, 1]))
  sdevs1 <- sqrt(t(x - m1) \%*\% s1inv \%*\% (x - m1))
  sdevs2 <- sqrt(t(x - m2) \%*\% s2inv \%*\% (x - m2))
  min(sdevs1, sdevs2)
}
implausibility <- function(x, target_level = 3) {
  sdtiny(x)
}
control_list <- list(
  box_limits = cbind(rep(-3, 2), rep(7, 2)),
  num_mutations = 40,
  num_iterations = 100
)
result <- subset_sim_slice_partitioned(
  implausibility,
  dims = 2,
  target_levels = 3,
  control_list = control_list,
  n_partitions = 2
)
}
