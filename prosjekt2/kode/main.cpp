#include <fstream>
#include <armadillo>
#include <cstdlib>
#include <cmath>
#include <ctime>

#include "jacobis_method.hpp"

int main(int argc, char * argv[]) {

  //deklarerer konstanter
  const int n = atoi(argv[1]); //leser inn dimensjonen n fra kommandolinja
  const double h = 1.0/n;
  const double hh = h*h;
  double a = -1/hh;
  double d = 2/hh;

  //deklarerer matriser og vektore
  mat A(n, n, fill::zeros);
  mat S(n, n, fill::eye);     //matrix to hold eigenvectors as row elements
  double *analytical_eigval = new double [n];
  vec jacobi_eigval(n);

  initialize(A, analytical_eigval, a, d, n);

  //bruker armadillo for å regne ut egenverdiene
  //tar tiden
  clock_t c_start = clock();
  vec arma_eigval = eig_sym(A); //gir egenverdiene i økende rekkefølge
  clock_t c_end = clock();

  // Beregner CPU-tid i milisekunder
  double arma_time = 1000.0 * (c_end - c_start) / CLOCKS_PER_SEC;

  //beregner egenvektorene med Jacobis metode
  //tar tiden
  int iterations;
  c_start = clock();
  jacobi(n, iterations, A, S, jacobi_eigval);
  c_end = clock();

  // Beregner CPU-tid i milisekunder
  double jacobi_time = 1000.0 * (c_end - c_start) / CLOCKS_PER_SEC;

  //kalkulerer maks relativ feil
  double max_error = calculate_max_error(n, arma_eigval, jacobi_eigval);

  write_data(n, iterations, arma_time, jacobi_time, max_error);

  return 0;
}
