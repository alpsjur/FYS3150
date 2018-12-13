#ifndef ROSSBY_PERIODIC_H
#define	ROSSBY_PERIODIC_H

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


void initWave(int, vec &, vector<double>&, bool);
void initialise(int, mat &, vec &, vector<double>);

void advance_vorticity_forward(double&, double, double, double, double);
void advance_vorticity_centered(double&, double, double, double, double, double);

void jacobisMethod2D(int, double, mat&, mat);

void writePsi(ofstream&, double&);
void writeZeta(ofstream&, double&);


inline int periodic(int j, int posdim){
  return (j + posdim) % posdim;
}
inline double sinewave(double x) {return sin(4.0*pi*x);}
inline double sinewaveDerivative(double x){return -16.0*pi*pi*sinewave(x);}
inline double gaussian(double x, double sigma) {return exp(-(x/sigma)*(x/sigma));}
inline double gaussianDerivative(double x, double sigma) {double sigma2 = sigma*sigma;
                    return -2.0*gaussian(x, sigma)*(sigma2 - 2.0*x*x)/(sigma2*sigma2);}

#endif /* ROSSBY_PERIODIC_H */
