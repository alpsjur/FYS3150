#include <iostream>

using namespace std;


void generate_elements(double& a[], double& b[], double& c[], double& d[]);
void solve_vector(double& v[],double& b[], double& c[], double& d[]);


int main(int argc, int * argv[]) {
  int n = argv[1];
  double * new v[n];
  double * new a[n];
  double * new b[n];
  double * new c[n];
  double * new d[n];

  generate_elements(double a, double b, double c, double d);
  delete[n] a;
  v[N] = d[N]/b[N];
  solve_vector(double v, double b, double c, double d);
  delete[n] b, c, d;
}


void generate_elements(double& a[], double& b[], double& c[], double& d[]) {
  for(i = 2; i < n; ++i) {
    double k = a[i-1]/b[i-1];
    b[i] -= k*c[i-1];
    d[i] -= k*d[i-1];
  }
}


void solve_vector(double& v[],double& b[], double& c[], double& d[]) {
  for(i = n-1; i > 1; --i) {
    v[i] = (d[i-1] - b[i-1]*v[i-1])/c[i-1];
  }
}
