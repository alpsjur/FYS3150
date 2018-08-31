#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>

using namespace std;


// deklarerer protofunksjoner som defineres senere
void init(int n, double*x, double *a, double *b, double *c, double *d);
void perform_gauss(int n, double *a, double *b, double *c, double *d);
void solve_vector(int n, double *v, double *b, double *c, double *d);
void calculate_error(int n, double *v, double *x, double *eps);
void write_data(int n, double *v, double *x, double *eps);

// definerer inline-funksjoner
inline double f(double xi) {return 100.0*exp(-10*xi);
}
inline double exact(double xi) {return 1.0 - (1.0 - exp(-10))*xi - exp(-10*xi);
}



int main(int argc, char * argv[]) {  // kommandolinje argumenter må være char
  const int n = atoi(argv[1]);     // std::atof : char -> int <cstdlib>

  // bruker pekere for dynamisk minne allokering, som slettes etter bruk
  try {
    // pekere for f, x, og feilen
    double *v = new double [n];
    double *x = new double [n];
    double *eps = new double [n];

    // pekere for matriseelementer
    double *a = new double [n];
    double *b = new double [n];
    double *c = new double [n];
    double *d = new double [n];

    init(n, x, a, b, c, d);
    perform_gauss(n, a, b, c, d);
    delete[] a;
    v[n-1] = d[n-1]/b[n-1];
    solve_vector(n, v, b, c, d);
    delete[] b;
    delete[] c;
    delete[] d;
    calculate_error(n, v, x, eps);
    write_data(n, v, x, eps);
    delete[] v;
  }
  catch(bad_alloc) {
    cout << "!FAILED TO ALLOCATE MEMORY FOR ARRAYS!" << '\n';
  }
  return 0;
}


// initialiserer matriseelementene og funksjonsverdiene
void init(int n, double *x, double *a, double *b, double *c, double *d) {
  const double h = 1.0/(n + 1.0);
  const double hh = h*h;

  for(int i = 0; i < n; ++i) {
    x[i] = (i + 1)*h;
    a[i] = -1.0;
    b[i] = 2.0;
    c[i] = -1.0;
    d[i] = hh*f(x[i]);  // eksponensialfunksjonen std::exp <cmath>
  }
}


// utfører gaussisk eliminasjon for å redusere likningssettet
void perform_gauss(int n, double *a, double *b, double *c, double *d) {
  for(int i = 1; i < n; ++i) {
    const double k = a[i-1]/b[i-1];
    b[i] -= k*c[i-1];
    d[i] -= k*d[i-1];
  }
}


// løser likningssettet for v-vektoren
void solve_vector(int n, double *v, double *b, double *c, double *d) {
  for(int i = n-2; i >= 0; --i) {
    v[i] = (d[i] - c[i]*v[i+1])/b[i];
  }
}


void calculate_error(int n, double *v, double *x, double *eps) {
  for(int i = 0; i < n; ++i) {
    eps[i] = log10(fabs((exact(x[i]) - v[i])/exact(x[i])));
  }
}

// skriver v-vektoren til .dat fil
// bruker pakken <fstream>
void write_data(int n, double *v, double *x, double *eps) {
  ofstream datafile;                // std::ofstream
  datafile.open("../data/general_matrix" + to_string(n) + ".dat");  // std::to_string
  for(int i = 0; i < n; ++i) {
    datafile << x[i] << ' ' <<v[i] << ' ' << eps[i] << endl;
  }
  datafile.close();
}
