
/*
løser likingssett på formen

  [b c 0 0 ... 0]   [v1]   [d1]
  [a b c 0 ... 0]   [v2]   [d2]
  [0 a b c ... 0]   [v3]   [d3]
  [0 0 a b ... 0] * [v4] = [d4]
  [......... b c]   [..]   [..]
  [0 ....... a b]   [vn]   [dn]

der a = c = -1, b = 2 og di = h^2*100*e^(-10xi)

*/
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <ctime>

using namespace std;


// deklarerer protofunksjoner som defineres senere
void init(int n, double*x, double *b, double *d);
void forward_sub(int n,  double *b, double *d);
void backward_sub(int n, double *v, double *b, double *d);
double calculate_error(int n, double *v, double *x, double *eps);
void write_data(int n, double *v, double *x, double *eps);
void write_error(int n, double max_error);
void write_CPU(int n, double CPU_time);

// definerer inline-funksjoner
inline double f(double xi) {return 100.0*exp(-10*xi);
}
inline double exact(double xi) {return 1.0 - (1.0 - exp(-10))*xi - exp(-10*xi);
}



int main(int argc, char * argv[]) {  // kommandolinje argumenter må være char
  const int n = pow(10, atoi(argv[1]));     // std::atof : char -> int <cstdlib>

  // bruker pekere for dynamisk minne allokering, som slettes etter bruk
  try {
    // pekere for f, x, og feilen
    double *v = new double [n];
    double *x = new double [n];
    double *eps = new double [n];

    // pekere for matriseelementer
    double *b = new double [n];
    double *d = new double [n];

    init(n, x, b, d);

    clock_t c_start = clock();
    forward_sub(n, b, d);
    v[n-1] = d[n-1]/b[n-1];
    backward_sub(n, v, b, d);
    clock_t c_end = clock();

    // Beregner CPU-tid i milisekunder
    double CPU_time = 1000.0 * (c_end - c_start) / CLOCKS_PER_SEC;
    write_CPU( n, CPU_time);

    delete[] b;
    delete[] d;
    double max_error = calculate_error(n, v, x, eps);

    if ( n < pow(10,4) ) {
      write_data(n, v, x, eps);
    }

    delete[] v;

    write_error( n, max_error);

  }
  catch(bad_alloc) {
    cout << "!FAILED TO ALLOCATE MEMORY FOR ARRAYS!" << '\n';
  }
  return 0;
}


// initialiserer matriseelementene og funksjonsverdiene
void init(int n, double *x, double *b, double *d) {
  const double h = 1.0/(n + 1.0);
  const double hh = h*h;

  b[0] = 2;
  x[0] = h;
  d[0] =  hh*f(x[0]);
  for (int i = 1; i<n; ++i){
    x[i] = (i+1)*h;
    b[i] = (double)(i+2)/(double)(i+1);
    d[i] = hh*f(x[i]);  // eksponensialfunksjonen std::exp <cmath>
  }
}


// utfører gaussisk eliminasjon for å redusere likningssettet
void forward_sub(int n,  double *b, double *d) {
  for(int i = 1; i < n; ++i) {
    d[i] += d[i-1]/b[i-1];
  }
}


// løser likningssettet for v-vektoren
void backward_sub(int n, double *v, double *b, double *d) {
  for(int i = n-2; i >= 0; --i) {
    v[i] = (d[i] + v[i+1])/b[i];
  }
}

double calculate_error(int n, double *v, double *x, double *eps) {
  double max_error = -10000.0;
  for(int i = 0; i < n; ++i) {
    eps[i] = log10(fabs((exact(x[i]) - v[i])/exact(x[i])));
    if (eps[i] > max_error){
      max_error = eps[i];
    }
  }
  return max_error;
}

// skriver v-vektoren til .dat fil
// bruker pakken <fstream>
void write_data(int n, double *v, double *x, double *eps) {
  ofstream datafile;                // std::ofstream
  datafile.open("../data/special_matrix" + to_string(n) + ".dat");  // std::to_string
  for(int i = 0; i < n; ++i) {
    datafile << x[i] << ' ' <<v[i] << ' ' << eps[i] << endl;
  }
  datafile.close();
}

// legger til n og max feil til .dat fil
void write_error(int n, double max_error) {
  ofstream logg;
  logg.open("../data/max_error_log.dat", fstream::app);
  logg << log10(n) << ' ' << max_error << endl;
  logg.close();
}

void write_CPU(int n, double CPU_time) {
  ofstream logg;
  logg.open("../data/special_matrix_time_log.dat", fstream::app);
  logg << log10(n) << ' ' << CPU_time << endl;
  logg.close();
}
