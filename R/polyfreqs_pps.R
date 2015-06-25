#' Posterior predictive model checks for polysomic inheritance
#'
#' @description Uses the posterior distribution of allele frequences from a \code{polyfreqs} run to test model fit using posterior predictive simulation.
#'
#' @param p_post A matrix containing the posterior samples from a \code{polyfreqs} run.
#' @param tM Total reads matrix: matrix containing the total number of reads mapping to each locus for each individual.
#' @param rM Reference reads marix: matrix containing the number of reference reads mapping to each locus for each individual.
#' @param ploidy Ploidy level of individuals in the population.
#' @param error The level of sequencing error. A fixed constant.
#'
#' @useDynLib polyfreqs
#' @importFrom Rcpp sourceCpp
#'
#' @export
polyfreqs_pps <- function(p_post, tM, rM, ploidy, error){

  sim_ref_read_ratios <- matrix(NA, nrow=nrow(p_post), ncol=ncol(p_post))
  sim_genos <- matrix(NA, nrow=nrow(tM),ncol=ncol(tM))
  sim_ref_read <- matrix(NA, nrow=nrow(tM),ncol=ncol(tM))
  obs_ref_read_ratio <- rep(NA, ncol(tM))
  obs_ref_read_ratio <- apply(rM/tM, 2, sum)

  for(i in 1:nrow(p_post)){
    sim_genos <- matrix(apply(as.matrix(p_post[i,]), 1, function(x) rbinom(nrow(tM), ploidy, x)), nrow=nrow(tM), ncol=ncol(tM))
    sim_ref_read <- sim_ref_reads(tM, sim_genos, ploidy, error)
    sim_ref_read_ratios[i,] <- apply(sim_ref_read/tM, 2, sum)
  }

  ratio_diff <- apply(sim_ref_read_ratios, 1, function(x) (obs_ref_read_ratio - x))
  locus_vec <- apply(ratio_diff, 2, function(x) sum(sign(quantile(x,c(0.025,0.975))))==0)
  return(list(ratio_diff=t(ratio_diff), locus_fit=locus_vec))
}
