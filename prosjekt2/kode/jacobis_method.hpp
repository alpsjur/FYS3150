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
void test_eigval(double *, vec, int);
double find_largest(mat, int *, int *, int);
void transform(mat &, mat &, int, int, int);
void write_data(int, int, double, double);

#endif /* JACOBIS_METHOD_H */
