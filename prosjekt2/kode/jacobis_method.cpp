#include <fstream>
#include <armadillo>
#include <cstdlib>
#include <cmath>
#include <ctime>

using namespace std;
using namespace arma;

void initialize(mat &A, double *analytical_eigval,double a, double d, int n);
double analytical(int i, int n, double a, double d);
void test_eigval(mat A, double *analytical_eigval, int n);
double find_largest(mat A, int *k, int *l, int n);
void transform(mat &A, mat &S,int k, int l, int n);

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

  test_eigval(A, analytical_eigval, n);

  double max = 10.0;
  double tol = 1e-10;
  int k, l, itterations;
  itterations = 0;
  while (max > tol){
    max = find_largest(A, &k, &l, n);
    transform(A, S, k, l, n);
    itterations++;
  }

  vec jacobi_eigval = A.diag(); //vektor med egenverdiene
  jacobi_eigval = sort(jacobi_eigval);
  cout << jacobi_eigval << endl;
  for (int i=0;i<n;i++){
    cout << analytical_eigval[n-1-i] << endl;
  }
  cout << itterations << endl;
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

void test_eigval(mat A, double *analytical_eigval, int n){
  /* Tester om egenverdiene til matrisen armadillo gir stemmer
  overens med de analytiske egenverdien
  */
  const double tol = 1e-3;

  //bruker armadillo for å regne ut egenverdiene
  vec arma_eigval = eig_sym(A); //gir egenverdiene i økende rekkefølge

  //itterer over egenverdiene og ser som differansen er innenfor tolleransen
  //arma_eigval er i stigende rekkefølge, analytical_eigval er i synkende
  for (int i=0; i<n;++i){
    if (fabs(arma_eigval(i)-analytical_eigval[n-1-i])>tol){
      cout << "Analytical and armadillo eigenvalues does not match" << endl;
      cout << "Diference is " << fabs(arma_eigval(i)-analytical_eigval[n-1-i]) << endl;
      break;
    }
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
        *k = i;  //ER DETTE RIKTIG?!?!?!
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