#include <Rcpp.h>
using namespace Rcpp;


// Sampling from full conditional for allele frequencies.

// [[Rcpp::export]]
NumericVector sample_p (IntegerMatrix Tm, IntegerMatrix Gn, SEXP pldy) {
  NumericVector pVnew(Tm.ncol());
  int ploidy = as<int>(pldy);
  for(int j = 0; j < Tm.ncol(); j++){
    checkUserInterrupt();
    int g_sum = 0;
    int gmin_sum = 0;
    for(int i = 0; i < Tm.nrow(); i++){
      if(Tm(i,j)==0){
        continue;
      }else{
        g_sum += Gn(i,j);
        gmin_sum += (ploidy-Gn(i,j));
      }
    }
    pVnew[j] = as<double>(rbeta(1,g_sum+1,gmin_sum+1));
  }
  return pVnew;
}
