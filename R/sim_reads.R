#' Simulation of sequencing read counts and genotypes
#'
#' @description Simulates genotypes and read counts under the model of Blischak et al. (2015).
#' @param pVec A vector of allele frequencies strung together using the concatenate function.
#' @param N_ind The number of individuals to simulate.
#' @param coverage The average number of sequences simulated per individual per locus (Poisson distributed).
#' @param ploidy The ploidy level of individuals in the population.
#' @param error The level of sequencing error. A fixed constant.
#'
#' @return A list of 3 matrices: (1) \code{genos} - a matrix of the simulated genotypes, (2) \code{tot_read_mat} - a matrix of the simulated number of total reads, (3) \code{ref_read_mat} - a matrix of the simulated number of reference reads.
#'
#' @references Blischak PD, Kubatko LS, Wolfe AD. 2015. Accounting for genotype uncertainty in the estimation of allele frequencies in autopolyploids. \emph{In review}. bioRxiv, \strong{doi}:####.
#' @useDynLib polyfreqs
#' @importFrom Rcpp sourceCpp

#' @export
sim_reads <- function(pVec, N_ind, coverage, ploidy, error){
  genos <- matrix(apply(as.matrix(pVec), 1, function(x) rbinom(N_ind, ploidy, x)), nrow=N_ind, ncol=length(pVec))
  tot_read_mat <- matrix(rpois(N_ind*length(pVec), coverage),nrow=N_ind, ncol=length(pVec))
  ref_read_mat <- sim_ref_reads(tot_read_mat, genos, ploidy, error)

  return(list(genos=genos, tot_read_mat=tot_read_mat, ref_read_mat=ref_read_mat))
}
