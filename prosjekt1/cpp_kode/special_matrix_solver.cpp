
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

using namespace std;

double error(double x, double sol);

int main(int argc, char *arg[])
{
  // trekker ut antall mesh points
  int n = atoi(arg[1]);

  // beregner step size
  double h = 1.0/(n+1);
  double hh = h*h;

  try
  {
    // initsialiserer d, b og x med dynamisk minnealkering
    double *d = new double[n];
    double *b = new double[n];
    double *x = new double[n];

    // beregner d og b
    d[0] = hh*100*exp(-10.0*h);
    b[0] = 2;
    x[0] = h;
    for (int i = 1; i<n; ++i)
    {
      x[i] = (i+1)*h;
      d[i] = hh*100*exp(-10.0*x[i]) + d[i-1]/b[i-1];
      b[i] = (double)(i+2)/(double)(i+1);
    }

    // initsialiserer sol med dynamisk minnealkoasjon
    double *v = new double[n];
    double *epsilon = new double[n];

    // beregner v og feilen
    v[n-1] = d[n-1]/b[n-1];
    epsilon[n-1] = error(x[n-1], v[n-1]);

    for (int i = n-2; i >= 0; --i)
    {
      v[i] = (d[i] + v[i+1])/b[i];
      epsilon[i] = error(x[i],sol[i]);
    }

    ofstream datafile;
    datafile.open("../data/special_matrix" + to_string(n) + ".dat"); 
    for(int i = 0; i < n; ++i)
    {
      datafile << x[i] << ' ' <<v[i] << ' ' << epsilon[i] << endl;
    }
    datafile.close();

    delete[] b;
    delete[] d;
    delete[] x;
    delete[] v;
    delete[] epsilon;
  }
    catch(bad_alloc) {
    cout << "!FAILED TO ALLOCATE MEMORY FOR ARRAYS!" << '\n';
  }
  return 0;
}

double error(double x, double sol)
{
  double exact = 1 - (1 - exp(-10)*x - exp(-10*x));
  double epsilon = log10(fabs(exact-sol)/exact);
  return epsilon;
}
