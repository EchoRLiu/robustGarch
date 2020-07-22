#include "RobGARCH_header_classes.h"
#include <cstdlib> // std::exit
//#include <numeric> 
#include <cmath> // std::sqrt(double), std::pow, std::abs

using std::cout;
using std::endl;
using std::exit;
using std::vector;
using std::array;
using std::sqrt;
using std::pow;
using std::accumulate;
using std::abs;

RobGarch11::RobGarch11(const vector<double>& x): x_(x)
{
  n_ = x.size();
  nDouble_ = n_ + 0.0;
  
  // get v_, vini_; normalize the data to be stored in y_.
  vX_ = tausq_(x);
  normalize_(); // Get normX_.
  vNX_ = tausq_(normX_);

}

void RobGarch11::sEst_(const vector<double> x)
{
  s_ = 1.0;
  double eps = 1.0;
  int nn = 1;
  
  Division xEmed(x, sEstEmed_);
  Division m(xEmed.abs());
  Division xDiv(x, m.mean());
  
}

vector<double> RobGarch11::rho_(const vector<double> x)
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
  double t = r.mean() * pow(s_, 2) / tausqConst_;
  return t;
}

void RobGarch11::normalize_(){
  
  double sqrtVX = sqrt(vX_);
  double meanNX = accumulate(x_.begin(), x_.end(), 0.0)/(nDouble_*sqrtVX);
  
  for(int i = 0; i < n_; ++i){
    normX_[i] = x_[i]/sqrtVX - meanNX; // Since there is an extra operation,
                                       // Division is not used.
    if(normX_[i] == 0){
      normX_[i] = zero_;
    }
  }
}

/////////////////////////////////////////////////////////
////////////////// Division Class ///////////////////////

// Constructors.
Division::Division(const vector<double> vecN): vecNumer_(vecN), denom_(1.0)
{
  n_ = vecN.size();
  nDouble_ = n_ + 0.0;
  getMean_();
  getAbs_();
}

Division::Division(const vector<double> vecN, const double d): vecNumer_(vecN), denom_(d)
{
  n_ = vecN.size();
  nDouble_ = n_ + 0.0;
  checkForZero_(d);
  getVecDiv_();
  getMean_();
  getAbs_();
}

// Accessors.
vector<double> Division::vecDiv() const
{
  return vecDiv_;
}

double Division::mean() const
{
  return mean_;
}

vector<double> Division::abs() const
{
  return vecDivAbs_;
}

vector<double> Division::vecNumer() const
{
  return vecNumer_;
}

double Division::denom() const
{
  return denom_;
}

// Mutators.

void Division::setVecNumer(const vector<double> vecNumer)
{
  vecNumer_ = vecNumer;
  
  n_ = vecNumer.size();
  nDouble_ = n_ + 0.0;

  getVecDiv_();
  getMean_();
  getAbs_();
}

void Division::setDenom(const double denom)
{
  checkForZero_(denom);
  denom_ = denom;
  getVecDiv_();
  getMean_();
  getAbs_();
}

// Private Member functions.

void Division::checkForZero_(double denom) const
{
  if(denom == 0.0)
  {
    std::exit(EXIT_FAILURE);
  }
}

void Division::getVecDiv_()
{
  for(int i = 0; i < n_; ++i)
  {
    vecDiv_[i] = vecNumer_[i] / denom_;
  }
}

void Division::getMean_()
{

  if(n_ == 0)
  { 
    std::exit(EXIT_FAILURE);
  }
  
  mean_ = 0.0;
  for(int i = 0; i < n_; ++i)
  {
    mean_ += vecNumer_[i];
  }
  mean_ /= nDouble_;
}

void Division::getAbs_()
{
  for(int i = 0; i < n_; ++i)
  {
    vecDiv_[i] = std::abs(vecDiv_[i]);
  }
}

Division Division::operator / (const double rhs) const
{
  Division div(vecDiv_, rhs);
  Division divDiv(div.vecDiv());
  
  return divDiv;
}

void Division::operator /= (const double rhs)
{
  *this = *this / rhs;
}

Division Division::operator - (const double rhs) const
{
  std::vector<double> vecDiff = vecDiv_;
  
  for(int i = 0; i < n_; i++)
  {
    vecDiff[i] -= rhs;
  }
  
  Division diff(vecDiff);
  
  return diff;
}

void Division::operator -= (const double rhs)
{
  *this = *this - rhs;
}
