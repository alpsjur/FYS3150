#ifndef ROSSBY_H
#define	ROSSBY_H

#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <ctime>

using namespace std;


void initMatrixElements(int&, double *, double *, double *, double *, double *);
void forward_sub(int&, double *, double *, double *, double *);
void backward_sub(int&, double *, double *, double *, double *);
void calculate_error(int&, double *, double *, double *);

void advance_vorticity_forward(double&, double&, double&, double&, double&, double&);
void advance_vorticity_centered(double&, double&, double&, double&, double&, double&);

void centered_difference(double&, double&, double&, double&);
void centered_difference2(double&, double&, double&, double&, double&);
void forward_difference(double&, double&, double&, double&);

inline double sinewave(double x) {return sin(4.0*pi*x);}
inline double gaussian(double x) {double sigma = 0.2; return exp((x/sigma)*(x/sigma));}

#endif /* ROSSBY_H */
