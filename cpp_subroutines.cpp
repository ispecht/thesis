#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
double epi(Rcpp::NumericVector q, double rho, double psi, bool obs, int N, int w, int M){
  double out = 0;
  if(q.length() > 0) {
    double lcoef = rho*log(psi) + log(1-psi) + log(rho);
    Rcpp::NumericVector qtmp (q.length(), 1.0);
    double lterm = lcoef;
    for(int j = 0; j < q.length(); ++j){
      qtmp[j] = qtmp[j] * (1 - q[j]);
      lterm = lterm + log(1 - qtmp[j]);
    }
    out = lterm;
    double ratio = 0;

    for(int i = 2; i <= M; ++i) {
      lcoef = lcoef + log(1-psi) + log(i + rho - 1) - log(i);
      lterm = lcoef;
      for(int j = 0; j < q.length(); ++j){
        qtmp[j] = qtmp[j] * (1 - q[j]);
        lterm = lterm + log(1 - qtmp[j]);
      }
      ratio = log(1 + (exp(lterm - out)));
      out = out + ratio;
    }
  }

  if(obs){
    out = out - log(N);
  }

  if(w > 0){
    out = epi(Rcpp::NumericVector (1, exp(out)), rho, psi, FALSE, N, w-1, M);
  }

  return out;
}



// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//
