#include "RobGARCH_header1.h"
#include <vector>
#include <cmath>
#include <algorithm>
#include <Rcpp.h>

using std::vector;
using std::copy;
using Rcpp::NumericVector;

// Tag every rcpp functions.

// [[Rcpp::export]]
NumericVector rcppbgarch11(NumericVector x){
  vector<double> stlVec(x.size());
  copy(x.begin(), x.end(), stlVec.begin());
  
  stlVec = bgarch11(stlVec);

  copy(stlVec.begin(), stlVec.end(), x.begin());
  return x;
}