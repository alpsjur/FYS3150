#ifndef ROSSBY_BOUNDED_H
#define	ROSSBY_BOUNDED_H

#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <iomanip>
#include <ctime>
#include <armadillo>

using namespace std;
using namespace arma;
double pi = 3.141592653589793;

void initWave(int, double, vector<double>&, vector<double>&, bool);

void initMatrixElements(int, double, vector<double>, vector<double>&, vector<double>&, vector<double>&, vector<double>&, vector<double>&);
void forward_sub(int, vector<double>&, vector<double>&, vector<double>&, vector<double>&);
void backward_sub(int, vector<double>&, vector<double>&, vector<double>&, vector<double>&);

void advance_vorticity_forward(double&, double, double, double, double);
void advance_vorticity_centered(double&, double, double, double, double, double);

void writePsi(ofstream&, double&);
void writeZeta(ofstream&, double&);

inline double sinewave(double x) {return sin(4.0*pi*x);}
inline double sinewaveDerivative(double x){return -16.0*pi*pi*sinewave(x);}
inline double gaussian(double x, double x0, double sigma) {return exp(-((x-x0)/sigma)*((x-x0)/sigma));}
inline double gaussianDerivative(double x, double x0,double sigma) {double sigma2 = sigma*sigma;
                    return -2.0*gaussian(x, x0, sigma)*(sigma2 - 2.0*(x-x0)*(x-x0))/(sigma2*sigma2);}

#endif /* ROSSBY_BOUNDED_H */
