#pragma
#ifndef ROBGARCH_HEADER_CLASSES_H
#define ROBGARCH_HEADER_CLASSES_H

#include <iostream>
#include <vector>
#include <array>

class RobGarch11
{
public:
  // Constructor: initialize x_, obtain vX_, normX_, vNX_.
  RobGarch11(const std::vector<double>& x);
  
  // prep for optimization in R - get v1, v2, vi, vini for Fnue.
  void prepOptBM1();
  // prep for optimization in R - get v1, v2, vi, vini for Fnue.
  void prepOptBM2();
  
  // Accessor for two vectors and two double for Fnue -> optimization in R.
  std::vector<double> getV1Fnue() const;
  std::vector<double> getV2Fnue() const;
  double getVi() const;
  double getVini() const;
  
  // The following is not needed anymore.
  // Final estimation is produced in R code.
  // Accessors: parameters estimations of BM1 and BM2.
  // std::array<double, 3> getBM1() const;
  // std::array<double, 3> getBM2() const;
  
private:
  
  void normalize_();
  double tausq_(const std::vector<double>& x);
  void sEst_(const std::vector<double> x);
  // future optimization: move rho_() to NewVector as its operation.
  std::vector<double> rho_(const std::vector<double> x) const;
  // More to come.
  
  
  std::vector<double> x_, normX_; // normX_ might contain const double zero_???
  double vX_, vNX_, s_;
  int n_;
  double nDouble_;
  
  const double zero_ = 10e-10;

  
  // More to come.
};



class NewVector
{
public:

  NewVector(const std::vector<double> vecN, const double d);
  NewVector(const std::vector<double> vecN);
  // Vector x is strongly preferred -> required?.
  // NewVector();
  
  // Access the NewVector result.
  std::vector<double> getVecDiv() const;
  double getMeanNumer() const; // mean of vecNumer_.
  std::vector<double> getAbs() const; //absolute value of vecDiv();
  std::vector<double> getVecNumer() const;
  double getDenom() const;
  
  // Mutators.
  void setVecNumer(const std::vector<double> vecNumer);
  void setDenom(const double denom);
  void setNumerDenom(const std::vector<double> vecNumer, const double denom);
  
  // Operators.
  NewVector operator / (const double rhs) const;
  void operator /= (const double rhs);
  NewVector operator - (const double rhs) const;
  void operator -= (const double rhs);
  
  // virtual ~NewVector() = default;

private:
  
  std::vector<double> vecNumer_, vecDiv_, vecDivAbs_; // ps_
  double denom_;
  double meanNumer_;
  int n_;
  double nDouble_;
  
  void checkForZero_(double d) const;
  void calcVecDiv_();
  void calcMeanNumer_(); // mean of vecNumer_.
  void calcAbs_();

};

class Extra
{
public:
  Extra(const std::vector<double> x);
  Extra(const std::vector<double> x, const double upper, const double lower);
  
  void setU(const double upper);
  void setL(const double lower);
  
  std::vector<double> getU() const;
  std::vector<double> getL() const;
  
  Extra operator * (std::vector<double> rhs) const;
  
private:
  
  std::vector<double> x_, u_, l_;
  double upperL_, lowerL_, nDouble_;
  int n_;
  
  void calcU_();
  void calcL_();
};

/*class RhoNewVector : public NewVector
{
  
};*/

#endif // ROBGARCH_HEADER_CLASSES_H

