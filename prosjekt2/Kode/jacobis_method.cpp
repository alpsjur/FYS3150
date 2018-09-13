#include <fstream>
#include <armadillo>
#include <cstdlib>
#include <cmath>
#include <ctime>

using namespace std;
using namespace arma;

void initialize(mat &A, double *analytical_eigval,double a, double d, int n);

//uttrykk for de analytiske egenverdiene
inline double analytical(int i, int n, double a, double d){
  return d+2*a*cos((i+1)*3.14159265/(n+1));
}

int main(int argc, char * argv[]) {

  //deklarerer konstanter
  const int n = atoi(argv[1]); //leser inn dimensjonen n fra kommandolinja
  const double h = 1.0/n;
  const double hh = h*h;
  double a = 2/hh;
  double d = -1/hh;

  //deklarerer matriser og vektorer
  mat A(n, n, fill::zeros);
  double *analytical_eigval = new double [n];

  initialize(A, analytical_eigval, a, d, n);

  //bruker armadillo for å regne ut egenverdiene
  vec arma_eigval = eig_sym(A); //gir egenverdiene i økende størrelse
  cout << arma_eigval << endl;
  cout << analytical_eigval << endl;

  return 0;
}

void initialize(mat &A, double *analytical_eigval,double a, double d, int n){
  //setter diagonalelementene til a og elementene over og under til d
  A(0,0) = a;
  analytical_eigval[0] = analytical(0, n, a, d);
  for (int i=1; i < n; ++i) {
    analytical_eigval[i] = analytical(i, n, a, d);
    A(i,i) = a;
    A(i-1,i) = d;
    A(i, i-1) = d;
  }
}
