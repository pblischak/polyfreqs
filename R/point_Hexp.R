#' Calculate expected heterozygosity from point estimates of allele frequencies
#'
#' @description Uses point estimates of allele frequencies calculated with \code{\link{simple_freqs}} to calculate the per locus and overall expected heterozygosity.
#' @param freqs A vector of estimated allele frequencies returned from the function \code{\link{simple_freqs}}.
#' @return A list with two items: the overall expected heterozygosity across loci (\code{Hexp}) and a vector with the per locus estimates of expected heterozygosity (\code{per_locus_Hexp}).

#' @export
point_Hexp <- function(freqs){
  allele1_sq <- freqs^2
  allele2_sq <- (1 - freqs)^2

  per_locus_sum <- allele1_sq + allele2_sq
  per_locus_Hexp <- 1 - per_locus_sum

  overall_sum <- sum(per_locus_sum)
  Hexp <- 1 - overall_sum/length(freqs)

  return(list(Hexp = Hexp, per_locus_Hexp = per_locus_Hexp))
}
