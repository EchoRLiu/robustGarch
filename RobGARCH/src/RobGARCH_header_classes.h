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
  
  // Parameters estimations of BM1.
  void optBM1();
  // Parameters estimations of BM2.
  void optBM2();
  
  // Accessors: parameters estimations of BM1 and BM2.
  std::array<double, 3> getBM1() const;
  std::array<double, 3> getBM2() const;
  
private:
  
  void normalize_();
  double tausq_(const std::vector<double>& x);
  void sEst_(const std::vector<double> x);
  std::vector<double> rho_(const std::vector<double> x);
  // More to come.
  
  
  std::vector<double> x_, normX_;
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
  
  NewVector operator / (const double rhs) const;
  void operator /= (const double rhs);
  NewVector operator - (const double rhs) const;
  void operator -= (const double rhs);

private:
  
  std::vector<double> vecNumer_, vecDiv_, vecDivAbs_;
  double denom_;
  double meanNumer_;
  int n_;
  double nDouble_;
  
  void checkForZero_(double d) const;
  void calcVecDiv_();
  void calcMeanNumer_(); // mean of vecNumer_.
  void calcAbs_();
  
};


#endif // ROBGARCH_HEADER_CLASSES_H

