#' Posterior predictive model checks for the model of inheritance
#'
#' @description Uses the posterior distribution of allele frequences from a \code{polyfreqs} run to test model fit using posterior predictive simulation.
#'
#' @param p_post A matrix containing the posterior samples from a \code{polyfreqs} run.
#' @param tM Total reads matrix.
#' @param rM Reference reads matrix.
#' @param ploidy Ploidy level of individuals in the population.
#' @param error Probability of a sequencing error.
#'
#' @useDynLib polyfreqs
#' @importFrom Rcpp sourceCpp
#'
#' @export
polyfreqs_pps <- function(p_post, tM, rM, ploidy, error){

  ref_read_liks <- matrix(NA,nrow=nrow(p_post),ncol=ncol(p_post))
  pred_genos <- matrix(apply(as.matrix(i), 1, function(x) rbinom(nrow(tM), ploidy, x)), nrow=nrow(tM), ncol=length(i))
  obs_ref_read_probs <- ref_read_pmf(tM, rM, pred_genos, error)
  obs_ref_read_lik <- apply(obs_ref_read_probs, 2, prod)

  for(i in 1:nrow(p_post)){
    pred_genos <- matrix(apply(as.matrix(i), 1, function(x) rbinom(nrow(tM), ploidy, x)), nrow=nrow(tM), ncol=length(i))
    pred_ref_reads <- sim_ref_reads(tM, pred_genos, ploidy, error)
    ref_read_probs <- ref_read_pmf(tM, pred_ref_reads, pred_genos, error)

    ref_read_liks[i,] <- apply(ref_read_probs, 2, prod)
  }
}
