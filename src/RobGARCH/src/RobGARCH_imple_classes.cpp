#include "RobGARCH_header_classes.h"
#include <cstdlib> // std::exit
//#include <numeric> 
#include <cmath> // std::sqrt(double), std::pow
#include <math.h> // std::fabs(double)

using std::cout;
using std::endl;
using std::exit;
using std::vector;
using std::array;
using std::sqrt;
using std::pow;
using std::accumulate;
using std::fabs;

RobGarch11::RobGarch11(const vector<double>& x): x_(x)
{
  n_ = x.size();
  nDouble_ = n_ + 0.0;
  
  // calculate v_, vini_; normalize the data to be stored in y_.
  vX_ = tausq_(x);
  normalize_(); // calculate normX_.
  vNX_ = tausq_(normX_);

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

double RobGarch11::tausq_(const vector<double>& x)
{
  const double tausqConst = 0.4797;
  
  sEst_(x);
  // In R code, s is set as global variable as "Sestim".
  // Here s_ is private member variable could be accessed.
  NewVector r(x, s_);
  r.setNumerDenom(rho_(r.getAbs()), 1.0);
  double t = r.getMeanNumer() * pow(s_, 2) / tausqConst;
  
  return t;
}

void RobGarch11::sEst_(const vector<double> x)
{
  const double sEstB = 1.625;
  const double sEstEmed = 0.675;
  const double sEstDiv1 = 0.405;
  const double sEstEpsComp = 0.005;
  
  s_ = 1.0;
  double eps = 1.0;
  int nn = 1; // n in R code.
  
  // m = median(fabs(x))/emed   //emed is positive.
  NewVector m(x, sEstEmed);
  m.setNumerDenom(m.getAbs(), 1.0); // m == m.getMeanNumer()
  // x /= m.
  NewVector xDiv(x, m.getMeanNumer()); // x = x/m == xDiv.getVecDiv()
  // rho1 = rho_(x/div1)
  // x/div1
  NewVector rho1(xDiv.getVecDiv(), sEstDiv1); //  x/div1 == rho1.getVecDiv()
  // rho_(x/div1)
  rho1.setNumerDenom(rho_(rho1.getAbs()), 1.0); // rho1 == rho1.getVecDiv()
  // a = mean(rho1)/b.
  double a = rho1.getMeanNumer() / sEstB;
  double v = 1-a;
  double si = a;
  // rho1 <- rho_(x/(sEstDiv1*si))
  // x / (sEstDiv1 * si)
  rho1 = xDiv / (sEstDiv1*si); // x/(sEstDiv1*si) == rho1.getVecDiv()
  // rho_(x/(sEstDiv1 * si))
  rho1.setNumerDenom(rho_(rho1.getAbs()), 1.0); // rho1 = rho_(x/(sEstDiv1*si)) == rho1.getVecDiv()
  a = rho1.getMeanNumer() / sEstB;
  double vi = 1-a;
  double AUX = v*vi;

  while(eps>sEstEpsComp && AUX>0.0)
  {
    nn += 1;
    s_ = si;
    v = vi;
    si = a * s_;
    rho1 = xDiv / (sEstDiv1*si);
    rho1.setNumerDenom(rho_(rho1.getAbs()), 1.0);
    a = rho1.getMeanNumer() / sEstB;
    vi = 1-a;
    AUX = v*vi;
    eps = std::fabs(s_ - si)/s_;
  }
  
  int nn2 = 0; // nsec in R code.
  double ns, nv;
  while(eps>sEstEpsComp)
  {
    ns = (s_ + si)/2.0;
    rho1 = xDiv / (sEstDiv1*ns);
    rho1.setNumerDenom(rho_(rho1.getAbs()), 1.0);
    a = rho1.getMeanNumer() / sEstB;
    nv = 1.0-a;
    AUX = nv * vi;
    
    if(AUX < 0.0)
    {
      v = nv; 
      s_ = ns;
    }
    else
    {
      vi = nv;
      si = ns;
    }
    eps = fabs(s_ - si) / s_;
    nn2 += 1;
    nn += 1;
  }
  
  s_ *= m.getMeanNumer();
  // if(nn > 30){ // n = n;}
  
}

vector<double> RobGarch11::rho_(const std::vector<double> x) const
{
  const double g1 = -1.944;
  const double g2 = 1.728;
  const double g3 = -0.312;
  const double g4 = 0.016;
  
  // Input x already has abs().
  Extra ax(x, 2.0, 3.0);
  
  
}


/////////////////////////////////////////////////////////
////////////////// NewVector Class ///////////////////////

// Constructors.
NewVector::NewVector(const vector<double> vecN): vecNumer_(vecN), denom_(1.0)
{
  n_ = vecN.size();
  nDouble_ = n_ + 0.0;
  vecDiv_ = vecN;
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

void NewVector::setNumerDenom(const vector<double> vecNumer, const double denom)
{
  checkForZero_(denom);
  vecNumer_ = vecNumer;
  denom_ = denom;
  
  n_ = vecNumer.size();
  nDouble_ = n_ + 0.0;
  
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
    vecDiv_[i] = fabs(vecDiv_[i]);
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



Extra::Extra(const vector<double> x): x_(x), upperL_(0.0), lowerL_(0.0)
{
  n_ = x.size();
  nDouble_ = n_ + 0.0;
  std::generate_n(std::back_inserter(u_), n_, [] { return 0.0; });
  std::generate_n(std::back_inserter(l_), n_, [] { return 0.0; });
  
  calcU_();
  calcL_();
}

void Extra::setU(const double upper)
{
  upperL_ = upper;
}

void Extra::setL(const double lower)
{
  lowerL_ = lower;
}

vector<double> Extra::getU() const
{
  return u_;
}
  
vector<double> Extra::getL() const
{
 return l_; 
}

Extra Extra::operator * (vector<double> rhs) const
{
  vector<double> prod(n_, 0.0);
  for(int i = 0; i < n_; i++)
  {
    prod[i] = x_[i] * rhs[i];
  }
  Extra prd(prod);
  
  return prd;
}

void Extra::calcL_()
{
  for(int i = 0; i < n_; i++)
  {
    if(x_[i] < upperL_){ l_[i] = 1.0;}
  }
}

void Extra::calcU_()
{
  for(int i = 0; i < n_; i++)
  {
    if(x_[i] > lowerL_){ u_[i] = 1.0;}
  }
}

/* vector<double> NewVector::operator * (const vector<double> rhs) const
 {
 vector<double> prod(n_, 0.0);
 
 for(int i = 0; i < n_; i++)
 {
 prod[i] = vecDiv_[i] * rhs[i];
 }
 
 return prod;
 } 

 std::vector<double> operator < (const double rhs) const;
 std::vector<double> operator > (const double rhs) const;
 
 
 vector<double> NewVector::operator < (const double rhs) const
 {
 vector<double> less(n_, 0.0); 
 
 for(int i = 0; i < n_; i++)
 {
 if(vecDiv_[i] < rhs){ less[i] = 1.0;}
 }
 
 return less;
 }
 
 
 vector<double> NewVector::operator > (const double rhs) const
 {
 vector<double> more(n_, 0.0);
 
 for(int i = 0; i < n_; i++)
 {
 if(vecDiv_[i] > rhs){ more[i] = 1.0;}
 }
 
 return more;
 }
 
 // r.setNumerDenom(rho_(r.getVecDiv()), 1.0);
 void NewVector::rho()
 {
 const double g1 = -1.944;
 const double g2 = 1.728;
 const double g3 = -0.312;
 const double g4 = 0.016;
 
 NewVector ax(vecDivAbs_);
 vector<double> u = ax > 3.0;
 vector<double> v = ax < 2.0;
 vector<double>
 
 }
 
 */





