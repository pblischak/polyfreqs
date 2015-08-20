#' Calculate posterior distribution of expected heterozygosity from posterior samples of allele frequencies
#'
#' @description Estimates a posterior distribution for the expected heterozygosity using the poterior samples of allele frequencies returned by \code{\link{polyfreqs}}.
#' @param p_post A matrix of posterior samples of allele frequencies returned by \code{\link{polyfreqs}}.
#' @param ploidy The ploidy level of individuals in the population (must be >= 2).
#' @return A list with two items: the posterior distribution for the overall expected heterozygosity across loci (\code{Hexp_post}) and a matrix with the per locus posterior samples of expected heterozygosity (\code{per_locus_Hexp_post}).

#' @export
posterior_Hexp <- function(p_post, ploidy){
  allele1_hom_mat <- p_post^ploidy
  allele2_hom_mat <- (1 - p_post)^ploidy

  per_locus_sum_mat <- allele1_hom_mat + allele2_hom_mat
  per_locus_Hexp_post <- 1 - per_locus_sum_mat

  overall_sum <- apply(per_locus_sum_mat, 1, sum)
  Hexp_post <- 1 - overall_sum/ncol(p_post)

  return(list(Hexp_post = Hexp_post, per_locus_Hexp_post = per_locus_Hexp_post))
}
