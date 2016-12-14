#include <Rcpp.h>

using namespace Rcpp;


int nonunif_int(const Rcpp::IntegerVector &vals, const Rcpp::NumericVector &probs){

  if(vals.size() != probs.size()){
    throw std::range_error("Vectors are not of the same length.");
  }

  double prob_sum = sum(probs);
  Rcpp::NumericVector normd_probs = probs/prob_sum;
  Rcpp::NumericVector cummulative_probs = Rcpp::cumsum(normd_probs);

  // Get random number between 0 and 1
  double samp = Rcpp::as<double>(Rcpp::runif(1,0,1));

  int res = -999;


  if(samp < cummulative_probs[0]){
    res = vals[0];
  } else {

    for(int i = 1; i < vals.size(); i++){

      if(samp >= cummulative_probs[i-1] && samp < cummulative_probs[i]){
        res = vals[i];
        break;
      }

    }
  }

  if(res == -999){
    Rcpp::stop("Invalid random integer generated");
  }

  return res;
}


// Full conditional sampling for genotypes

// [[Rcpp::export]]
IntegerMatrix sample_g(IntegerMatrix Tm, IntegerMatrix Rm, IntegerMatrix Gn, NumericVector pV, SEXP pldy, SEXP err) {

  int ploidy = as<int>(pldy);
  double error = as<double>(err);
  double gVec_sum = 0.0, gEpsilon = 0.0;
  IntegerMatrix GnNew(Tm.nrow(),Tm.ncol());
  NumericVector gVec_tmp(ploidy+1);
  NumericVector gVec(ploidy+1);
  IntegerVector gCount(ploidy+1);

  for(int r = 0; r <= ploidy; r++){
    gCount[r] = r;
  }

  for(int j = 0; j < Tm.ncol(); j++){
    checkUserInterrupt();
    for(int i = 0; i< Tm.nrow(); i++){
      // Check for missing data in total read matrix (i.e. entry equal 0)
      if(Tm(i,j)==0){
        continue;
      }else{
        gVec_sum = 0.0;
        for(int p = 0; p <= ploidy; p++){


          gEpsilon = p / (double) ploidy * (1 - error)
                   + (1 - p / (double) ploidy) * error;

          gVec_tmp[p] = ::Rf_dbinom(Rm(i,j), Tm(i,j), gEpsilon, 0)
                      * ::Rf_dbinom(p, ploidy, pV[j], 0);

          gVec_sum += gVec_tmp[p];
          //if(p==0){
          //  gVec_tmp[p] = pow(error,Rm(i,j))*pow((1-error),(Tm(i,j)-Rm(i,j)))*pow((1-pV[j]),ploidy);
          //}else if(p==ploidy){
          //  gVec_tmp[p] = pow((1-error),Rm(i,j))*pow((error),(Tm(i,j)-Rm(i,j)))*pow(pV[j],ploidy);
          //}else{
          //  ratio = p/(double)ploidy;
          //  ratio_w_error = ratio*(1-error) + (1-ratio)*error;
          //  gVec_tmp[p] = Rf_choose(ploidy,p)*pow(ratio_w_error,Rm(i,j))*pow((1-ratio_w_error),(Tm(i,j)-Rm(i,j)))*pow(pV[j],p)*pow((1-pV[j]),(ploidy-p));
          //}

        }


        //gVec_sum = std::accumulate(gVec_tmp.begin(),gVec_tmp.end(),0.0);

        for(int d = 0; d <= ploidy; d++){
          gVec[d] = gVec_tmp[d] / gVec_sum;
          /*Rcout << d << "\t" << gVec_tmp[d] << "\t"
                  << gVec[d] << "\t" << gVec_sum << std::endl;*/
        }

        GnNew(i,j) = nonunif_int(gCount, gVec);

      }
    }
  }

  return GnNew;

}
