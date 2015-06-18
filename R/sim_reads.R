#' @useDynLib polyfreqs
#' @importFrom Rcpp sourceCpp

#' @export
sim_reads <- function(pVec, N_ind, coverage, ploidy, error){
  genos <- matrix(apply(as.matrix(pVec), 1, function(x) rbinom(N_ind, ploidy, x)), nrow=N_ind, ncol=length(pVec))
  tot_read_mat <- matrix(rpois(N_ind*length(pVec), coverage),nrow=N_ind, ncol=length(pVec))
  ref_read_mat <- sim_ref_reads(tot_read_mat, genos, ploidy, error)

  return(list(genos=genos, tot_read_mat=tot_read_mat, ref_read_mat=ref_read_mat))
}
