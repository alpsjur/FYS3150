
/*
løser likingssett på formen

  [b c 0 0 ... 0]   [v0]   [d0]
  [a b c 0 ... 0]   [v1]   [d1]
  [0 a b c ... 0]   [v2]   [d2]
  [0 0 a b ... 0] * [v3] = [d3]
  [......... b c]   [..]   [..]
  [0 ....... a b]   [vn]   [dn]

der a = c = -1, b = 2 og di = h^2*100*e^(-10i*h) der h = 1/n

*/

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>

using namespace std;

void write_data(int n, double *v);

int main(int argc, char *arg[])
{
  // trekker ut antall mesh points
  int n = atoi(arg[1]);

  // beregner step size
  const double h = 1.0/(n+1);
  const double hh = h*h;

  try
  {
    // initsialiserer d og b med dynamisk minnealkoasjon
    double *d = new double[n];
    double *b = new double[n];

    // beregner d og b
    d[0] = hh*100*exp(-10.0*h);
    b[0] = 2;
    for (int i = 1; i<n; ++i)
    {
      d[i] = hh*100*exp(-10.0*(i+1)*h) + d[i-1]/b[i-1];
      b[i] = (double)(i+2)/(double)(i+1);
    }

    // initsialiserer v med dynamisk minnealkoasjon
    double *v = new double[n];

    // beregner v
    v[n-1] = d[n-1]/b[n-1];
    for (int i = n-2; i >= 0; --i)
    {
      v[i] = (d[i] + v[i+1])/b[i];
    }
    write_data(n, v);
    delete[] b;
    delete[] d;
    delete[] v;
  }
    catch(bad_alloc) {
    cout << "!FAILED TO ALLOCATE MEMORY FOR ARRAYS!" << '\n';
  }
  return 0;
}

// skriver v-vektoren til .dat fil
void write_data(int n, double *v)
{
  ofstream datafile;
  datafile.open("../data/special_matrix" + to_string(n) + ".dat");  // std::to_string
  for(int i = 0; i < n; ++i)
  {
    datafile << v[i] << endl;
  }
  datafile.close();
}
