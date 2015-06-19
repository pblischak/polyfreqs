#include <Rcpp.h>
using namespace Rcpp;

// Probability mass function for the number of observed reference reads.
// Also the likelihood (Eq. 2 in the manuscript).

// [[Rcpp::export]]
NumericMatrix ref_read_pmf(IntegerMatrix Tm, IntegerMatrix Rm, IntegerMatrix Gn, SEXP pldy, SEXP err){

  double error = as<double>(err);
  int ploidy = as<int>(pldy);
  double geno_ratio;
  NumericMatrix RefReadProbs(Tm.nrow(),Tm.ncol());

  for(int i = 0; i < Tm.nrow(); i++){
    for(int j = 0; j < Tm.ncol(); j++){

      // Check for missing data in total read matrix (i.e. entry equal 0)
      if(Tm(i,j)==0){
        continue;
      }else{
        if(Gn(i,j)==0){
          RefReadProbs(i,j) = R::dbinom(Rm(i,j), Tm(i,j), error, FALSE);
        }else if(Gn(i,j)==ploidy){
          RefReadProbs(i,j) = R::dbinom(Rm(i,j), Tm(i,j), 1-error, FALSE);
        }else{
          geno_ratio = Gn(i,j)/(double)ploidy;
          RefReadProbs(i,j) = R::dbinom(Rm(i,j), Tm(i,j), geno_ratio, FALSE);
        }
      }
    }
  }
  return RefReadProbs;
}

