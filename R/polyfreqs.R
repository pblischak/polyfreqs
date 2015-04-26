#' Bayesian estimation of allele frequencies in polyploids
#'
#' @description \code{polyfreqs} implements a Gibbs sampler to perform Bayesian inference on the allele frequencies in a population of polyploids.
#' @author Paul Blischak
#' @param tM Total reads matrix: matrix containing the total number of reads mapping to each locus for each individual.
#' @param rM Reference reads marix: matrix containing the number of reference reads mapping to each locus for each individual.
#' @param ploidy The ploidy level of the sampled individuals (must be >= 2).
#' @param iter The number of MCMC generations to run (default=50,000).
#' @param thin Thins the MCMC output by sampling everything \code{thin} generations (default=100).
#' @param error The level of sequencing error. A fixed constant (default=0.01).
#' @param genotypes Logical variable indicating whether or not to print the values of the genotypes sampled during the MCMC (default=FALSE).
#' @param col_header Optional column header tag for use in running loci in parallel (default="").
#' @param outfile The name of the ouput file that samples from the posterior distribution of allele frequencies are written to (default="polyfreqs-mcmc.out").
#' @examples
#' data(total_reads)
#' data(ref_reads)
#' polyfreqs(total_reads,ref_reads,4,iter=100,thin=10)
#'
#' @useDynLib polyfreqs
#' @importFrom Rcpp sourceCpp
#' @import RcppArmadillo

#' @export
polyfreqs<- function(tM, rM, ploidy, iter=50000, thin=100, error=0.01, genotypes=FALSE, col_header="", outfile="polyfreqs-mcmc.out"){

  # Check that input matrices are valid.
  stopifnot(is.matrix(tM))
  stopifnot(is.matrix(rM))

  if(genotypes){
    d.success<-dir.create("./genotypes")
    if(!d.success){
      stop("Attempt to make genotypes directory failed.")
    }

    # Print column headers for the output files.
    cat("iter", paste("p_",col_header,"_loc",1:ncol(tM),sep=""), sep="\t", file=outfile)
    cat("\n",file=outfile, append=TRUE)

    ## Initialize genotype matrix from discrete uniform(0,ploidy)
    ## Initialize allele frequency vector from uniform(0,1)
    ## Replace entries in genotype matrix with missing data in
    ## total read and reference read files (tM = 0).
    missing.data<-(tM==0)
    gM_init<-matrix(sample(0:ploidy,nrow(tM)*ncol(tM),replace=TRUE), nrow(tM), ncol(tM))
    gM_init[missing.data]=0
    pV_init<-runif(ncol(tM))



    rnames<-row.names(tM)
    for(i in 1:nrow(tM)){
      cat("iter", paste("g_",col_header,"_loc",1:ncol(tM),sep=""),sep="\t", file=paste("./genotypes/",rnames[i],"_g-mcmc.out",sep=""))
      cat("\n",file=paste("./genotypes/", rnames[i],"_g-mcmc.out",sep=""), append=TRUE)
    }

    # Start MCMC
    for(k in 1:iter){

      # Sample from full conditional on genotypes first.
      # Then reassign the values.
      gM<-sample_g(tM, rM, gM_init, pV_init, ploidy, error)
      gM_init<-gM

      # Sample from the full conditional on allele frequencies
      # given draw from genotypes.
      pV<-sample_p(tM, gM_init, ploidy)
      pV_init<-pV

      # Print every 'thin' generation of the MCMC.
      if(k %% thin == 0){
        cat(k, pV_init, sep="\t",file=outfile, append=TRUE)
        cat("\n",file=outfile, append=TRUE)
        print_g(k,gM_init,tM)
      }
    }

  } else {
    # Print column headers for the output files.
    cat("iter", paste("p_",col_header,"_loc",1:ncol(tM),sep=""), sep="\t", file=outfile)
    cat("\n",file=outfile, append=TRUE)

    ## Initialize genotype matrix from discrete uniform(0,ploidy)
    ## Initialize allele frequency vector from uniform(0,1)
    ## Replace entries in genotype matrix with missing data in
    ## total read and reference read files (tM = 0).
    missing.data<-(tM==0)
    gM_init<-matrix(sample(0:ploidy,nrow(tM)*ncol(tM),replace=TRUE), nrow(tM), ncol(tM))
    gM_init[missing.data]=0
    pV_init<-runif(ncol(tM))


    # Start MCMC
    for(k in 1:iter){

      # Sample from full conditional on genotypes first.
      # Then reassign the values.
      gM<-sample_g(tM, rM, gM_init, pV_init, ploidy, error)
      gM_init<-gM

      # Sample from the full conditional on allele frequencies
      # given draw from genotypes.
      pV<-sample_p(tM, gM_init, ploidy)
      pV_init<-pV

      # Print every 'thin' generation of the MCMC.
      if(k %% thin == 0){
        cat(k, pV_init, sep="\t",file=outfile, append=TRUE)
        cat("\n",file=outfile, append=TRUE)
      }
    }

  }

}
