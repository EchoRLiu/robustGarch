#include "RobGARCH_header_classes.h"
#include <cstdlib> // std::exit
//#include <numeric> 
#include <cmath> // std::sqrt(double), std::pow

using std::cout;
using std::endl;
using std::exit;
using std::vector;
using std::array;
using std::sqrt;
using std::pow;
using std::accumulate;

RobGarch11::RobGarch11(const vector<double>& x): x_(x)
{
  
  // get v_, vini_; normalize the data to be stored in y_.
  vX_ = tausq_(x);
  normalize_(); // Get normX_.
  vNX_ = tausq_(normX_);

}

void RobGarch11::sEst_(const vector<double>& x)
{
  
  // More to come.
}

vector<double> RobGarch11::rho_(const vector<double>& x)
{
  
  // More to come.
}

double RobGarch11::tausq_(const vector<double>& x)
{
  sEst_(x);
  // In R code, s is set as global variable as "Sestim".
  // Here s_ is private member variable could be accessed.
  Division xS(x, s_);
  Division r(rho_(xS.vecDiv()));
  double t = r.mean() * std::pow(s_, 2) / tausqConst_;

}

void RobGarch11::normalize_(){
  
  double sqrtVX = sqrt(vX_);
  double meanNX = accumulate(x_.begin(), x_.end(), 0.0)/(x_.size()*sqrtVX);
  
  for(int i = 0; i < x_.size(); ++i){
    normX_[i] = x_[i]/sqrtVX - meanNX; // Since there is an extra operation,
                                       // Division is not used.
    if(normX_[i] == 0){
      normX_[i] = 10.0e-10;
    }
  }
}

/////////////////////////////////////////////////////////
////////////////// Division Class ///////////////////////

Division::Division(const vector<double> vecN): vecNumer_(vecN), denom_(1.0)
{
  getMean_();
}

Division::Division(const vector<double> vecN, const double d): vecNumer_(vecN), denom_(d)
{
  checkForZero_(d);
  getVecDiv_();
  getMean_();
}

vector<double> Division::vecDiv() const
{
  return vecDiv_;
}

double Division::mean() const
{
  return mean_;
}

vector<double> Division::vecNumer() const
{
  return vecNumer_;
}

double Division::denom() const
{
  return denom_;
}

void Division::setVecNumer(const vector<double> vecNumer)
{
  vecNumer_ = vecNumer;
  getVecDiv_();
}

void Division::setDenom(const double denom)
{
  checkForZero_(denom);
  denom_ = denom;
  getVecDiv_();
}

void Division::checkForZero_(double denom) const
{
  if(denom == 0.0)
  {
    std::exit(EXIT_FAILURE);
  }
}

void Division::getVecDiv_()
{
  for(int i = 0; i < vecNumer_.size(); ++i)
  {
    vecDiv_[i] = vecNumer_[i] / denom_;
  }
}

void Division::getMean_()
{
  double n = vecNumer_.size() + 0.0;
  if(n == 0.0)
  { 
    std::exit(EXIT_FAILURE);
  }
  
  mean_ = 0.0;
  for(int i = 0; i < vecNumer_.size(); ++i)
  {
    mean_ += vecNumer_[i];
  }
  mean_ /= n;
}


