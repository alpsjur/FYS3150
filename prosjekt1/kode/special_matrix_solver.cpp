
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

double error(double x, double approx);

int main(int argc, char *arg[])
{
  // trekker ut antall mesh points
  int n = pow(10,atoi(arg[1]));

  // beregner step size
  double h = 1.0/(n+1);
  double hh = h*h;

  try
  {
    // initsialiserer d, b og x med dynamisk minnealokering
    double *b = new double[n];
    double *d = new double[n];
    double *x = new double[n];

    b[0] = 2;
    x[0] = h;
    d[0] =  hh*100*exp(-10.0*x[0]);
    for (int i = 1; i<n; ++i)
    {
      x[i] = (i+1)*h;
      b[i] = (double)(i+2)/(double)(i+1);
      d[i] = hh*100.0*exp(-10.0*x[i]);
    }

    // initsialiserer v og d med dynamisk minnealkoasjon
    double *v = new double[n];

    clock_t c_start = clock();
    // forward substitution
    for (int i = 1; i < n; ++i)
    {
      d[i] =  d[i] - d[i-1]/b[i-1];
    }

    // backward substitution
    v[n-1] = d[n-1]/b[n-1];
    for (int i = n-2; i >= 0; --i)
    {
      v[i] = (d[i] + v[i+1])/b[i];
    }
    clock_t c_end = clock();

    cout << 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC << endl;

    delete[] b;
    delete[] d;

    /*
    ofstream datafile;
    datafile.open("../data/special_matrix" + to_string(n) + ".dat");
    for(int i = 0; i < n; ++i)
    {
      datafile << x[i] << ' ' <<v[i] << ' ' << error(x[i], v[i]) << endl;
    }
    datafile.close()
    */
    delete[] x;
    delete[] v;
  }
    catch(bad_alloc) {
    cout << "!FAILED TO ALLOCATE MEMORY FOR ARRAYS!" << '\n';
  }
  return 0;
}

double error(double x, double approx)
{
  double exact = 1.0 - (1.0 - exp(-10.0))*x - exp(-10.0*x);
  double epsilon = log10(fabs((exact-approx)/exact));
  return epsilon;
}
