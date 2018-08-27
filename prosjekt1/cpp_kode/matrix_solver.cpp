#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;


void initialise_arrays(int n, double& f[], double h);
void generate_gauss(int n,double& a[], double& b[], double& c[], double& d[]);
void solve_vector(int n, double& v[],double& b[], double& c[], double& d[]);
void write_data(int n, double v[]);


int main(int argc, int * argv[]) {
  const int n = argv[1];
  const int h = 1/(n + 1);

  double * new v[n];
  double * new a[n];
  double * new b[n];
  double * new c[n];
  double * new d[n];

  initialise_arrays(n, b, h);
  generate_gauss(n, a, b, c, d);
  delete[n] a;
  v[n] = d[n]/b[n];
  solve_vector(n, v, b, c, d);
  delete[n] b, c, d;

  return 0;
}


void initialise_arrays(int n, double& b[], double h) {
  for(int i = 0; i < n+1; ++i) {
    a[i] = -1.0;
    b[i] = 2.0;
    c[i] = -1.0;
    d[i] = h**2*100*exp(-10*i*h);

  }
}


void generate_gauss(int n, double& a[], double& b[], double& c[], double& d[]) {
  for(int i = 2; i < n; ++i) {
    double k = a[i-1]/b[i-1];
    b[i] -= k*c[i-1];
    d[i] -= k*d[i-1];
  }
}


void solve_vector(int n, double& v[],double& b[], double& c[], double& d[]) {
  for(int i = n-1; i > 1; --i) {
    v[i] = (d[i-1] - b[i-1]*v[i-1])/c[i-1];
  }
}


void write_data(int n, double v[]) {
  ofstream datafile;
  datafile.open("../data/matrix" << n << ".dat");
  datafile << v;
  datafile.close();
}
