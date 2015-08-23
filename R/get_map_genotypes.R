#' Maximum \emph{a posteriori} (MAP) estmation of autopolyploid genotypes
#'
#' @description Calculates the MAP estimate of the genotypes for autopolyploid individuals.
#' @param tM Total reads matrix: matrix containing the total number of reads mapping to each locus for each individual.
#' @param burnin Percent of the posterior samples to discard as burn-in (default=20).
#' @param geno_dir File path to directory containing the posterior samples of genotypes output by \code{\link{polyfreqs}} (default = "genotypes").
#' @return A matrix containing the maximum \emph{a posteriori} estimates for all individuals at each locus. The MAP estimate of the genotype is simply the posterior mode.

#' @export
get_map_genotypes <- function(tM, burnin=20, geno_dir="genotypes"){
  names <- rownames(tM)
  nind <- nrow(tM)
  nloci <- ncol(tM)
  if(!(length(names)==nind)){
    stop("The number of names and individuals does not match.")
  }

  map_genotypes <- matrix(NA, nrow = nind, ncol = nloci)

  Mode <- function(g) {
    g_unique <- unique(g)
    g_unique[which.max(tabulate(match(g, g_unique)))]
  }

  for(i in 1:length(names)){
    tmp_table <- read.table(paste("./",geno_dir,"/",names[i],"_g-mcmc.out",sep=""), header=T, row.names=1)
    tmp_vec <- apply(as.matrix(tmp_table[round(burnin/100 * nrow(tmp_table) + 1):nrow(tmp_table),]), 2, Mode)
    map_genotypes[i,] <- tmp_vec
  }

  return(map_genotypes)
}
