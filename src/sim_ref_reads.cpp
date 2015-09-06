#include <Rcpp.h>
using namespace Rcpp;

// Function for simulation from the likelihood function
// for reference sequencing reads. Used for simulating data
// and for posterior predictive simulations.

// [[Rcpp::export]]
IntegerMatrix sim_ref_reads(IntegerMatrix Tm, IntegerMatrix Gn, SEXP pldy, SEXP err) {

  int ploidy = as<int>(pldy);
  double error = as<double>(err);
  double geno_ratio;
  double geno_ratio_error;
  IntegerMatrix RefMat(Tm.nrow(),Tm.ncol());

  for(int i = 0; i < Tm.nrow(); i++){
    for(int j = 0; j < Tm.ncol(); j++){
      // Check for missing data in total read matrix (i.e. entry equal 0)
      if(Tm(i,j)==0){
        continue;
      }else{
        if(Gn(i,j)==0){
          RefMat(i,j) = as<int>(rbinom(1, Tm(i,j), error));
        }else if(Gn(i,j)==ploidy){
          RefMat(i,j) = as<int>(rbinom(1, Tm(i,j), 1-error));
        }else{
          geno_ratio = (double)Gn(i,j)/(double)ploidy;
          geno_ratio_error = geno_ratio*(1-error) + (1-geno_ratio)*error;
          RefMat(i,j) = as<int>(rbinom(1, Tm(i,j), geno_ratio_error));
        }
      }
    }
  }
  return RefMat;
}
