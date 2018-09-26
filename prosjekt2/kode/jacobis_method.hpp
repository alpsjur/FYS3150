#ifndef JACOBIS_METHOD_H
#define	JACOBIS_METHOD_H

#include <fstream>
#include <armadillo>
#include <cstdlib>
#include <cmath>
#include <ctime>

using namespace std;
using namespace arma;

void initialize(mat &, double *, double, double, int);
double analytical(int, int, double, double);
double find_largest(mat, int &, int &, int);
void transform(mat &, mat &, int, int, int);
void jacobi(int, int &, mat, mat &, vec &);
void write_data(int, int, double, double, double);
double calculate_max_error(int, vec, vec);

inline double analytical(int i, int n, double a, double d){
  const double pi = 3.14159;
  return d+2.0*a*cos((i+1)*pi/(n+1.0));
}

#endif /* JACOBIS_METHOD_H */
