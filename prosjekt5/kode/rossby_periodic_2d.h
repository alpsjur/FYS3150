#ifndef ROSSBY_PERIODIC_2D_H
#define	ROSSBY_PERIODIC_2D_H

#include <fstream>
#include <armadillo>
#include <cstdlib>
#include <iomanip>
#include <cmath>
#include <ctime>
#include <vector>

using namespace std;
using namespace arma;

double pi = 3.1415926535897932;

void advance_vorticity_forward(double&, double, double, double, double);
void advance_vorticity_centered(double&, double, double, double, double, double);

void writePsi(ofstream&, double&);
void writeZeta(ofstream&, double&);


inline int periodic(int j, int posdim){
  return (j + posdim) % posdim;
}
inline double sinewave(double x, double y) {return sin(4.0*pi*x+4.0*pi*x);}
inline double sinewaveDerivative(double x, double y){return -16.0*pi*pi*(sinewave(x)+sinewave(y));}
inline double gaussian(double x, double y, double sigma) {return exp(-(1/2/sigma)*(1/2/sigma)*(x*x+y*y));}
inline double gaussianDerivative(double x, double sigma) {double sigma2 = sigma*sigma;
                    return -2.0*gaussian(x, sigma)*(sigma2 - 2.0*x*x)/(sigma2*sigma2);}

#endif /* ROSSBY_PERIODIC_2D_H */
