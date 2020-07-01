// Header file. Should include all functions in reusable cpp files.
#include <vector>
using std::vector;

vector<double> bgarch11(vector<double> x);

double varin(vector<double> x);

double tausq(vector<double> x);

double SEST2(vector<double> x);

vector<double> RHO2(vector<double> x);

vector<double> nnest1(vector<double> y, double vini);

double Fnue(vector<double> vi, double vini);

vector<double> nfun(vector<double> x);

vector<double> freg(vector<double> x, double a, double b);

vector<double> gk(vector<double> x, double k, double l);

vector<double> nnest2(vector<double> y, double vini);

double Fnue2(vector<double> vi, double vini);

vector<double> nfun2(vector<double> x);