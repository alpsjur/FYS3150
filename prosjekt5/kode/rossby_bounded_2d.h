#ifndef ROSSBY_BOUNDED_2D_H
#define	ROSSBY_BOUNDED_2D_H

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

void initWave(int, double, mat&, mat&, bool);

void advance_vorticity_forward(double&, double, double, double, double);
void advance_vorticity_centered(double&, double, double, double, double, double);

void jacobisMethod2D(int, double, mat&, mat);

void writePsi(ofstream&, double&);
void writeZeta(ofstream&, double&);


inline double sinewave(double x, double y) {return sin(4.0*pi*x+4.0*pi*x);}
inline double sinewaveDerivative(double x, double y){return -16.0*pi*pi*2*sinewave(x,y);}
inline double gaussian(double x, double x0, double y, double y0, double sigma) {return exp(-(1/sigma)*(1/sigma)/2*((x-x0)*(x-x0)+(y-y0)*(y-y0)));}
inline double gaussianDerivative(double x, double x0, double y, double y0, double sigma) {double sigma2 = sigma*sigma*2;
                    return (gaussian(x, x0, y, y0, sigma)/(sigma*sigma/2)/(sigma2*sigma2))*4*(x*x+y*y-2*sigma);}

#endif /* ROSSBY_BOUNDED_2D_H */
