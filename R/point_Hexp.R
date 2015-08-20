#' Calculate expected heterozygosity from point estimates of allele frequencies
#'
#' @description Uses point estimates of allele frequencies calculated with \code{\link{simple_freqs}} to calculate the per locus and overall expected heterozygosity.
#' @param freqs A vector of estimated allele frequencies returned from the function \code{\link{simple_freqs}}.
#' @param ploidy The ploidy level of individuals in the population (must be >= 2).
#' @return A list with two items: the overall expected heterozygosity across loci (\code{Hexp}) and a vector with the per locus estimates of expected heterozygosity (\code{per_locus_Hexp}).

#' @export
point_Hexp <- function(freqs, ploidy){
  allele1_hom <- freqs^ploidy
  allele2_hom <- (1 - freqs)^ploidy

  per_locus_sum <- allele1_hom + allele2_hom
  per_locus_Hexp <- 1 - per_locus_sum

  overall_sum <- sum(per_locus_sum)
  Hexp <- 1 - overall_sum/length(freqs)

  return(list(Hexp = Hexp, per_locus_Hexp = per_locus_Hexp))
}
