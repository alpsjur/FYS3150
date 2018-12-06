#ifndef ROSSBY_H
#define	ROSSBY_H

#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <iomanip>
#include <ctime>
#include <armadillo>

using namespace std;
double pi = 3.1416;

void initPsi(int, vector<double>&);

void initMatrixElements(int, vector<double>, vector<double>&, vector<double>&, vector<double>&, vector<double>&, vector<double>&);
void forward_sub(int, vector<double>&, vector<double>&, vector<double>&, vector<double>&);
void backward_sub(int, vector<double>&, vector<double>&, vector<double>&, vector<double>&);


void advance_vorticity_forward(double&, double&, double&, double&, double&, double&);
void advance_vorticity_centered(double&, double&, double&, double&, double&, double&);

void writePsi(ofstream&, double&);
void writeZeta(ofstream&, double&);

void centered_difference(double&, double&, double&, double&);
void centered_difference2(double&, double&, double&, double&, double&);
void forward_difference(double&, double&, double&, double&);


inline double sinewave(double x) {return sin(4.0*pi*x);}
inline double gaussian(double x) {double sigma = 0.2; return exp((x/sigma)*(x/sigma));}

#endif /* ROSSBY_H */