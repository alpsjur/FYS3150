#include <fstream>
#include <armadillo>
#include <cstdlib>
#include <cmath>
#include <ctime>

#include "jacobis_method.hpp"

int main(int argc, char * argv[]) {

  //deklarerer konstanter
  //leser inn dimensjonen n fra kommandolinja
  const int n = atoi(argv[1]);

  // 0 = buckling beam, 1 = 1 quantum dot, 2 = 2 quantum dots
  const int problem = atoi(argv[2]);

  // 0 = skriver ikke egenparene til fil, 1 = skriver egenparene til fil
  const int eigwrite = atoi(argv[3]);

  // a  trengs hvis vi vil teste eigenverdiene med de analytiske
  double a;
  int iterations;

  // tar omega_r fra terminalen hvis vi ser på to elektroner
  double omega_r = 1;
  int interact = 0;
  if(problem == 2){
    omega_r = atof(argv[4]);
    interact = atof(argv[5]);
  }

  double rhomin = 0.0;
  double rhomax = 10.0;

  double *d = new double [n];

  //deklarerer matriser og vektorer
  mat A(n, n, fill::zeros);
  mat S(n, n, fill::eye);     //matrix to hold eigenvectors as row elements
  vec jacobi_eigval(n);

  initialize(A, d, a, rhomin, rhomax, omega_r, problem, interact, n);

  //bruker armadillo for å regne ut egenverdiene og egenvektorene
  vec arma_eigval;
  mat arma_eigvec;
  //tar tiden
  clock_t c_start = clock();
  eig_sym(arma_eigval, arma_eigvec, A);  //gir egenverdiene i økende rekkefølge
  clock_t c_end = clock();

  // Beregner CPU-tid i milisekunder
  double arma_time = 1000.0 * (c_end - c_start) / CLOCKS_PER_SEC;

  //beregner egenverdiene og egenvektorene med Jacobis metode
  //tar tiden
  c_start = clock();
  jacobi(n, iterations, A, S, jacobi_eigval);
  c_end = clock();

  // Beregner CPU-tid i milisekunder
  double jacobi_time = 1000.0 * (c_end - c_start) / CLOCKS_PER_SEC;

  //kalkulerer maks relativ feil
  double max_error = calculate_max_error(n, a, d, problem,jacobi_eigval);

  delete[] d;

  write_log(n, iterations, arma_time, jacobi_time, max_error);
  if(eigwrite == 1){
    write_eig(n,problem, jacobi_eigval, S);
  }

  return 0;
}
