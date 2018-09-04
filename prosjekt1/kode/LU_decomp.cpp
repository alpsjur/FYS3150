#include <fstream>
#include <armadillo>
#include <cstdlib>
#include <cmath>
#include <ctime>

using namespace std;
using namespace arma;


void initialise(int n, mat A, colvec f);
void test_lu(int n, mat L, mat U, mat P, mat A);
void write_data(int n, colvec x);

inline double init_f(double xi) {return 100.0*exp(-10*xi);
}


int main(argc, char *argv[]) {
  const int n = pow(10, atoi(argv[1]));

  // declaring matrices
  mat A = <mat>(n, n).zeros(); // laplacian matrix for Dirichlet bc
  mat L = <mat>(n, n);         // lower triangular matrix
  mat U = <mat>(n, n);         // upper triangular matrix
  mat P = <mat>(n, n);         // permutation matrix

  // declaring vectors
  colvec f = <colvec>(n);      // column vector containing function values
  colvec x;                    // column vector in eq Ux = y
  colvec y;                    // column vector in eq Ly = f

  initialise(n, A, f);         // initialising A with tridiagonal values
  lu(L, U, P, A);              // performing LU-decomposition on A
  test_lu(n , L, U, P, A);     // checking to see that matrix mult returns A
  solve(y, trimatu(L), f);     // solving for y indicating that L is triangular
  solve(x, trimatu(U), y);     // solving for x, our solution
  write_data(n, x);

  return 0;
}


void initialise(int n, mat A, colvec f) {
  const double h = 1.0/(n + 1);
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


void test_lu(int n, mat L, mat U, mat P, mat A) {
  mat B = P.t()*L*U;
  mat diff = A - B;
  rowsum = sum(diff, 1);
  colsum = sum(diff);
  if(fabs(rowsum) > 1e-6 and fabs(colsum) > 1e-6) {
    cout << "!A - LU != 0!" << endl;
  }
}

void write_data(int n, colvec x) {
  ofstream datafile;                // std::ofstream
  datafile.open("../data/LU_decomp" + to_string(n) + ".dat");  // std::to_string
  for(int i = 0; i < n; ++i) {
    datafile << x[i] << endl;
  }
  datafile.close();
}
