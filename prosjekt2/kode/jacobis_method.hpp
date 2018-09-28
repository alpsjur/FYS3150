#ifndef JACOBIS_METHOD_H
#define	JACOBIS_METHOD_H

#include <fstream>
#include <armadillo>
#include <cstdlib>
#include <cmath>
#include <ctime>

using namespace std;
using namespace arma;

void initialize(mat &, double *, double &, double, double, double, int, int);
double analytical(int, int, double, double);
double find_largest(mat, int &, int &, int);
void transform(mat &, mat &, int, int, int);
void jacobi(int, int &, mat, mat &, vec &);
void write_log(int, int, double, double, double);
double calculate_max_error(int, vec, double *);
void write_eig(int, int,vec, mat);

inline double analytical_buck(int i, int n, double a, double di){
  const double pi = 3.14159;
  return di+2.0*a*cos((i+1)*pi/(n+1.0));
}

inline double analytical_dot(int i){
  return 3+i*4;
}

#endif /* JACOBIS_METHOD_H */
