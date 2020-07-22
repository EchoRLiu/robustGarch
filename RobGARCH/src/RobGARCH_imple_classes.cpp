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
  
  // calculate v_, vini_; normalize the data to be stored in y_.
  vX_ = tausq_(x);
  normalize_(); // calculate normX_.
  vNX_ = tausq_(normX_);

}

void RobGarch11::sEst_(const vector<double> x)
{
  const double sEstB = 1.625;
  const double sEstEmed = 0.675;
  
  s_ = 1.0;
  double eps = 1.0;
  int nn = 1;
  
  NewVector xEmed(x, sEstEmed);
  NewVector m(xEmed.getAbs());
  NewVector xDiv(x, m.getMeanNumer());
  
}

vector<double> RobGarch11::rho_(const vector<double> x)
{
  
  // More to come.
}

double RobGarch11::tausq_(const vector<double>& x)
{
  const double tausqConst = 0.4797;
  
  sEst_(x);
  // In R code, s is set as global variable as "Sestim".
  // Here s_ is private member variable could be accessed.
  NewVector xS(x, s_);
  NewVector r(rho_(xS.getVecDiv()));
  double t = r.getMeanNumer() * pow(s_, 2) / tausqConst;
  return t;
}


// Done.
void RobGarch11::normalize_(){
  
  double sqrtVX = sqrt(vX_);
  double meanNX = accumulate(x_.begin(), x_.end(), 0.0)/(nDouble_*sqrtVX);
  
  for(int i = 0; i < n_; ++i){
    normX_[i] = x_[i]/sqrtVX - meanNX; // Since there is an extra operation,
                                       // NewVector is not used.
    if(normX_[i] == 0){
      normX_[i] = zero_;
    }
  }
}

/////////////////////////////////////////////////////////
////////////////// NewVector Class ///////////////////////

// Constructors.
NewVector::NewVector(const vector<double> vecN): vecNumer_(vecN), denom_(1.0)
{
  n_ = vecN.size();
  nDouble_ = n_ + 0.0;
  calcMeanNumer_();
  calcAbs_();
}

NewVector::NewVector(const vector<double> vecN, const double d): vecNumer_(vecN), denom_(d)
{
  n_ = vecN.size();
  nDouble_ = n_ + 0.0;
  checkForZero_(d);
  calcVecDiv_();
  calcMeanNumer_();
  calcAbs_();
}

// Accessors.
vector<double> NewVector::getVecDiv() const
{
  return vecDiv_;
}

double NewVector::getMeanNumer() const
{
  return meanNumer_;
}

vector<double> NewVector::getAbs() const
{
  return vecDivAbs_;
}

vector<double> NewVector::getVecNumer() const
{
  return vecNumer_;
}

double NewVector::getDenom() const
{
  return denom_;
}

// Mutators.

void NewVector::setVecNumer(const vector<double> vecNumer)
{
  vecNumer_ = vecNumer;
  
  n_ = vecNumer.size();
  nDouble_ = n_ + 0.0;

  calcVecDiv_();
  calcMeanNumer_();
  calcAbs_();
}

void NewVector::setDenom(const double denom)
{
  checkForZero_(denom);
  denom_ = denom;
  calcVecDiv_();
  calcMeanNumer_();
  calcAbs_();
}

// Private Member functions.

void NewVector::checkForZero_(double denom) const
{
  if(denom == 0.0)
  {
    std::exit(EXIT_FAILURE);
  }
}

void NewVector::calcVecDiv_()
{
  for(int i = 0; i < n_; ++i)
  {
    vecDiv_[i] = vecNumer_[i] / denom_;
  }
}

void NewVector::calcMeanNumer_()
{

  if(n_ == 0)
  { 
    std::exit(EXIT_FAILURE);
  }
  
  meanNumer_ = 0.0;
  for(int i = 0; i < n_; ++i)
  {
    meanNumer_ += vecNumer_[i];
  }
  meanNumer_ /= nDouble_;
}

void NewVector::calcAbs_()
{
  for(int i = 0; i < n_; ++i)
  {
    vecDiv_[i] = std::abs(vecDiv_[i]);
  }
}

NewVector NewVector::operator / (const double rhs) const
{
  NewVector div(vecDiv_, rhs);
  NewVector divDiv(div.getVecDiv());
  
  return divDiv;
}

void NewVector::operator /= (const double rhs)
{
  *this = *this / rhs;
}

NewVector NewVector::operator - (const double rhs) const
{
  std::vector<double> vecDiff = vecDiv_;
  
  for(int i = 0; i < n_; i++)
  {
    vecDiff[i] -= rhs;
  }
  
  NewVector diff(vecDiff);
  
  return diff;
}

void NewVector::operator -= (const double rhs)
{
  *this = *this - rhs;
}
