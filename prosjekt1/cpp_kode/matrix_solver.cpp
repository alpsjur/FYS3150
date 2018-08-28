#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>

using namespace std;


// deklarerer protofunksjoner som defineres senere
void init(int n, double *a, double *b, double *c, double *d, double h);
void perform_gauss(int n, double *a, double *b, double *c, double *d);
void solve_vector(int n, double *v, double *b, double *c, double *d);
void write_data(int n, double *v);



int main(int argc, char * argv[]) {  // kommandolinje argumenter må være char
  const int n = atof(argv[1]);     // std::atof : char -> int <cstdlib>
  const double h = 1.0/(n + 1.0);

  // bruker pekere for dynamisk minne allokering, som slettes etter bruk
  try {
    double *v = new double [n];
    double *a = new double [n];
    double *b = new double [n];
    double *c = new double [n];
    double *d = new double [n];

    init(n, a, b, c, d, h);
    perform_gauss(n, a, b, c, d);
    delete[] a;
    v[n-1] = d[n-1]/b[n-1];
    solve_vector(n, v, b, c, d);
    delete[] b;
    delete[] c;
    delete[] d;
    write_data(n, v);
    delete[] v;
  }
  catch(bad_alloc) {
    cout << "failed to allocate memory!" << '\n';
  }
  return 0;
}


// initialiserer matriseelementene og funksjonsverdiene
void init(int n, double *a, double *b, double *c, double *d, double h) {
  for(int i = 1; i < n; ++i) {
    a[i] = -1.0;
    b[i] = 2.0;
    c[i] = -1.0;
    d[i] = h*h*100.0*exp(-10.0*i*h);  // eksponensialfunksjonen std::exp <cmath>
  }
}


// utfører gaussisk eliminasjon for å redusere likningssettet
void perform_gauss(int n, double *a, double *b, double *c, double *d) {
  for(int i = 2; i < n; ++i) {
    const double k = a[i-1]/b[i-1];
    b[i] -= k*c[i-1];
    d[i] -= k*d[i-1];
  }
}


// løser likningssettet for v-vektoren
void solve_vector(int n, double *v, double *b, double *c, double *d) {
  for(int i = n-2; i > 0; --i) {
    v[i] = (d[i] - c[i]*v[i+1])/b[i];
  }
}


// skriver v-vektoren til .dat fil
// bruker pakken <fstream>
void write_data(int n, double *v) {
  ofstream datafile;
  datafile.open("../data/matrix" + to_string(n) + ".dat");  // std::to_string
  for(int i = 0; i < n+1; ++i) {
    datafile << v[i] << endl;
  }
  datafile.close();
}
