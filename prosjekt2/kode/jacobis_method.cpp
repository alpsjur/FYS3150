#include <fstream>
#include <armadillo>
#include <cstdlib>
#include <cmath>
#include <ctime>

using namespace std;
using namespace arma;

void initialize(mat &A, double *analytical_eigval,double a, double d, int n);
double analytical(int i, int n, double a, double d);
void test_eigval(double *analytical_eigval, vec arma_eigval, int n);
double find_largest(mat A, int *k, int *l, int n);
void transform(mat &A, mat &S,int k, int l, int n);
void write_data(int n, int itterations, double arma_time, double jacobi_time);

int main(int argc, char * argv[]) {

  //deklarerer konstanter
  const int n = atoi(argv[1]); //leser inn dimensjonen n fra kommandolinja
  const double h = 1.0/n;
  const double hh = h*h;
  double a = 2/hh;
  double d = -1/hh;

  //deklarerer matriser og vektorer
  mat A(n, n, fill::zeros);
  mat S(n, n, fill::eye);     //matrix to hold eigenvectors as row elements
  double *analytical_eigval = new double [n];

  initialize(A, analytical_eigval, a, d, n);

  //bruker armadillo for å regne ut egenverdiene
  //tar tiden
  clock_t c_start = clock();
  vec arma_eigval = eig_sym(A); //gir egenverdiene i økende rekkefølge
  clock_t c_end = clock();

  // Beregner CPU-tid i milisekunder
  double arma_time = 1000.0 * (c_end - c_start) / CLOCKS_PER_SEC;

  test_eigval(analytical_eigval, arma_eigval, n);

  //beregner egenvektorene med Jacobis metode
  double max = 10.0;
  double tol = 1e-10;
  int k, l, itterations;
  itterations = 0;
  c_start = clock();   //tar tiden
  while (max > tol){
    max = find_largest(A, &k, &l, n);
    transform(A, S, k, l, n);
    itterations++;
  }
  c_end = clock();

  // Beregner CPU-tid i milisekunder
  double jacobi_time = 1000.0 * (c_end - c_start) / CLOCKS_PER_SEC;

  vec jacobi_eigval = A.diag(); //vektor med egenverdiene
  jacobi_eigval = sort(jacobi_eigval);

  write_data(n, itterations, arma_time, jacobi_time);

  return 0;
}

void initialize(mat &A, double *analytical_eigval,double a, double d, int n){
  //setter diagonalelementene til a og elementene direkte over og under til d
  A(0,0) = d;
  analytical_eigval[0] = analytical(0, n, a, d);
  for (int i=1; i < n; ++i) {
    analytical_eigval[i] = analytical(i, n, a, d);
    A(i,i) = d;
    A(i-1,i) = a;
    A(i, i-1) = a;
  }
  return;
}

double analytical(int i, int n, double a, double d){
  const double pi = 3.14159;
  return d+2.0*a*cos((i+1)*pi/(n+1.0));
}

void test_eigval(double *analytical_eigval, vec arma_eigval, int n){
  /*
  Tester om egenverdiene til matrisen armadillo gir stemmer
  overens med de analytiske egenverdien
  */
  const double tol = 1e-3;
  int count = 0;
  double tot_rel_err = 0;

  //itterer over egenverdiene og ser som relativ feil er innenfor tolleransen
  //arma_eigval er i stigende rekkefølge, analytical_eigval er i synkende
  for (int i=0; i<n;++i){
    if (fabs((arma_eigval[i]-analytical_eigval[n-1-i])/
        analytical_eigval[n-1-i])>tol){
      count++;
      tot_rel_err += fabs((arma_eigval[i]-analytical_eigval[n-1-i])/
                     analytical_eigval[n-1-i]);
    }
  }
  if (count > 0){
    cout << "Analytisk og armadillo egenverdier var forskellige " <<
            count << " ganger for n=" << n << endl;
    cout << "I snitt var relativ feil " << tot_rel_err/count << endl;
  }
  return;
}

double find_largest(mat A, int *k, int *l, int n) {
  double largest = 0.0;
  //ittererer over øvre halvdel av matrisen, siden matrisen er symmetrisk
  for (int i=0; i<n; ++i){
    for (int j=i+1;j<n;++j){
      double aij = fabs(A(i,j));
      if (aij>largest){
        largest = aij;
        *k = i;
        *l = j;
      }
    }
  }
  return largest;
}

void transform(mat &A, mat &S,int k, int l, int n) {
  double tau, t;
  tau = (A(l,l) - A(k,k))/(2.0*A(k,l));

  //for å unngå loss of precission når t -> inf skriver vi om uttrykket
  if (tau >= 0){
    t = 1.0/(tau+sqrt(1.0+tau*tau));   //tangens theta
  }
  else {
    t = -1.0/(-tau+sqrt(1.0+tau*tau)); //tangens theta
  }

  double c, cc, s, ss, cs;
  c = 1.0/(sqrt(1+t*t));             //cosinus theta
  cc = c*c;
  s = c*t;                           //sinus theta
  ss = s*s;
  cs = c*s;

  //oppdaterer A
  double akk;
  akk = A(k,k);
  A(k,k) = akk*cc - 2*A(k,l)*cs + A(l,l)*ss;
  A(l,l) = A(l,l)*cc + 2*A(k,l)*cs + akk*ss;
  A(k,l) = 0.0;//(akk-all)*cs + akl*(cc-ss);
  A(l,k) = 0.0;//A(k,l);
  for (int i=0; i<n; ++i) {
    if (i != k && i != l) {
      double aik;
      aik = A(i,k);
      A(i,k) = aik*c - A(i,l)*s;
      A(k,i) = A(i,k);
      A(i,l) = A(i,l)*c + aik*s;
      A(l,i) = A(i,l);
    }
    //kalkulerer egenvektorene
    double sik, sil;
    sik = S(i,k);
    sil = S(i,l);
    S(i,k) = c*sik - s*sil;
    S(i,l) = c*sil + s*sik;
  }
  return;
}

void write_data(int n, int itterations, double arma_time, double jacobi_time){
  ofstream logg;
  logg.open("../data/jacobi_log.dat", fstream::app);
  logg << n << ' ' << itterations << ' ' << jacobi_time << ' ' << arma_time << endl;
  logg.close();
}
