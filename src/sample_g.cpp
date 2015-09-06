// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadilloExtensions/sample.h>

using namespace Rcpp;


// [[Rcpp::export]]
IntegerVector csample_integer( IntegerVector x, int size, bool replace,
NumericVector prob = NumericVector::create()) {
  RNGScope scope;
  IntegerVector ret = RcppArmadillo::sample(x, size, replace, prob);
  return ret;
}


// Full conditional sampling for genotypes

// [[Rcpp::export]]
IntegerMatrix sample_g(IntegerMatrix Tm, IntegerMatrix Rm, IntegerMatrix Gn, NumericVector pV, SEXP pldy, SEXP err) {

  int ploidy = as<int>(pldy);
  int num_samps = 1;
  double error = as<double>(err);
  bool rep = TRUE;
  double ratio;
  double ratio_w_error;
  IntegerMatrix GnNew(Tm.nrow(),Tm.ncol());
  NumericVector gVec_tmp(ploidy+1);
  NumericVector gVec(ploidy+1);
  IntegerVector gCount(ploidy+1);

  for(int r = 0; r <= ploidy; r++){
    gCount[r] = r;
  }

  for(int j = 0; j < Tm.ncol(); j++){
    for(int i = 0; i< Tm.nrow(); i++){
      // Check for missing data in total read matrix (i.e. entry equal 0)
      if(Tm(i,j)==0){
        continue;
      }else{
        for(int p = 0; p <= ploidy; p++){
          if(p==0){
            gVec_tmp[p] = pow(error,Rm(i,j))*pow((1-error),(Tm(i,j)-Rm(i,j)))*pow((1-pV[j]),ploidy);
          }else if(p==ploidy){
            gVec_tmp[p] = pow((1-error),Rm(i,j))*pow((error),(Tm(i,j)-Rm(i,j)))*pow(pV[j],ploidy);
          }else{
            ratio = p/(double)ploidy;
            ratio_w_error = ratio*(1-error) + (1-ratio)*error;
            gVec_tmp[p] = Rf_choose(ploidy,p)*pow(ratio_w_error,Rm(i,j))*pow((1-ratio_w_error),(Tm(i,j)-Rm(i,j)))*pow(pV[j],p)*pow((1-pV[j]),(ploidy-p));
          }

        }

        double gVec_sum = 0;
        gVec_sum = std::accumulate(gVec_tmp.begin(),gVec_tmp.end(),0.0);

        for(int d = 0; d <= ploidy; d++){
          gVec[d] = gVec_tmp[d]/gVec_sum;
        }

        GnNew(i,j) = as<int>(csample_integer(gCount,num_samps,rep,gVec));
      }
    }
  }
  return GnNew;
}
