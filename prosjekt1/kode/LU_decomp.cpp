#include <fstream>
#include <armadillo>
#include <cstdlib>
#include <cmath>
#include <ctime>

using namespace std;
using namespace arma;


void initialise(int n, mat &A, vec &f);
void write_data(int n, vec x, double CPU_time);

inline double init_f(double xi) {
  return 100.0*exp(-10*xi);
}


int main(int argc, char *argv[]) {
  const int n = pow(10, atoi(argv[1]));

  // declaring matrices
  mat A(n, n, fill::zeros); // laplacian matrix for Dirichlet bc
  mat L;         // lower triangular matrix
  mat U;         // upper triangular matrix

  // declaring vectors
  vec f(n);      // column vector containing function values
  vec y(n);         // solution of Ly = f
  vec x(n);         // solution of Ux = y

  initialise(n, A, f);         // initialising A with tridiagonal values
  clock_t c_start = clock();
  lu(L, U, A);              // performing LU-decomposition on A
  solve(y, trimatl(L), f);     // solving for y indicating that L is triangular
  solve(x, trimatu(U), y);     // solving for x, our solution
  clock_t c_end = clock();

  // Beregner CPU-tid i milisekunder
  double CPU_time = 1000.0 * (c_end - c_start) / CLOCKS_PER_SEC;

  write_data(n, x, CPU_time);

  return 0;
}


void initialise(int n, mat &A, vec &f) {
  const double h = 1.0/(n + 1.0);
  const double hh = h*h;
  f[0] = hh*init_f(0);
  A(0,0) = 2.0;
  for (int i=1; i < n; ++i) {
    f[i] = hh*init_f(i*h);
    A(i,i) = 2.0;
    A(i-1,i) = -1.0;
    A(i, i-1) = -1.0;
  }
}


void write_data(int n, vec x, double CPU_time) {
  ofstream datafile;                // std::ofstream
  datafile.open("../data/LU_decomp" + to_string(n) + ".dat");  // std::to_string
  for(int i = 0; i < n; ++i) {
    datafile << x[i] << endl;
  }
  datafile.close();

  ofstream logg;
  logg.open("../data/LU_time_log.dat", fstream::app);
  logg << log10(n) << ' ' << CPU_time << endl;
  logg.close();
}
