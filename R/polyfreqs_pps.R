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

  #sim_ref_read_liks <- matrix(NA, nrow=nrow(p_post), ncol=ncol(p_post))
  sim_ref_read_ratios <- matrix(NA, nrow=nrow(p_post), ncol=ncol(p_post))
  #obs_ref_read_liks <- matrix(NA,nrow=nrow(p_post),ncol=ncol(p_post))
  obs_ref_read_ratio <- apply(rM/tM, 2, sum)

  for(i in 1:nrow(p_post)){
    sim_genos <- matrix(apply(p_post[i,], 2, function(x) rbinom(nrow(tM), ploidy, x)), nrow=nrow(tM), ncol=length(p_post[i,]))
    sim_ref_read <- sim_ref_reads(tM, sim_genos, ploidy, error)
    sim_ref_read_ratios[i,] <- apply(sim_ref_read/tM, 2, sum)
    sim_ref_read_probs <- ref_read_pmf(tM, sim_ref_read, sim_genos, ploidy, error)
    obs_ref_read_probs <- ref_read_pmf(tM, rM, sim_genos, ploidy, error)
    #sim_ref_read_liks[i,] <- apply(log(sim_ref_read_probs), 2, sum)
    #obs_ref_read_liks[i,] <- apply(log(obs_ref_read_probs), 2, sum)
  }

  #likelihood_diff <- obs_ref_read_liks - sim_ref_read_liks
  ratio_diff <- apply(sim_ref_read_ratios, 1, function(x) (obs_ref_read_ratio - x))
  locus_vec <- apply(ratio_diff, 2, function(x) sum(sign(quantile(x,c(0.025,0.975))))==0)
  return(list(ratio_diff=ratio_diff, locus_fit=locus_vec))
}
