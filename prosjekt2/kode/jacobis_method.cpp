#include <fstream>
#include <armadillo>
#include <cstdlib>
#include <cmath>
#include <ctime>

#include "jacobis_method.hpp"


void jacobi(int n, int &iterations, mat A, mat &S, vec &jacobi_eigval){
  double max = 10.0;
  double tol = 1e-10;
  int k, l;
  iterations = 0;
  while (max > tol){
    max = find_largest(A, k, l, n);
    transform(A, S, k, l, n);
    iterations++;
  }
  jacobi_eigval = A.diag();            //vektor med egenverdiene
  return;
}


void initialize(mat &A, double *d, double &a, double rhomin, double rhomax, double omega_r, int problem, int interact, int n){
  // definerer step size og diagonalelementer
  const double h = (rhomax - rhomin)/double(n);
  const double hh = h*h;
  a = -1.0/hh;

  for(int i = 0; i < n; ++i){
    d[i] = 2.0/hh;

    double rho = (rhomin + (i+1)*h);
    double V = rho*rho;
    // sjekker hvilket problem vi ser på, og velger potensialet
    if(problem == 1){
      d[i] += V;
    }
    else if (problem == 2){
      d[i] += omega_r*V;
      if(interact == 1){
        d[i] += 1.0/rho;
      }
    }
  }
  //setter diagonalelementene i A til d og elementene direkte over og under til a
  A(0,0) = d[0];
  for (int i=1; i < n; ++i) {
    A(i,i) = d[i];
    A(i-1,i) = a;
    A(i, i-1) = a;
  }
  return;
}

double find_largest(mat A, int &k, int &l, int n) {
  //finner indeksene til det største elementet i matrisen, og returnerer maxverdien
  double largest = 0.0;
  //ittererer over øvre halvdel av matrisen, siden matrisen er symmetrisk
  for (int i=0; i<n; ++i){
    for (int j=i+1; j<n; ++j){
      double aij = fabs(A(i,j));
      if (aij>largest){
        largest = aij;
        k = i;
        l = j;
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
  c = 1.0/(sqrt(1+t*t));               //cosinus theta
  cc = c*c;
  s = c*t;                             //sinus theta
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
  //oppdaterer S
    double sik;
    sik = S(i,k);
    S(i,k) = c*sik - s*S(i,l);
    S(i,l) = c*S(i,l) + s*sik;
  }
  return;
}


double calculate_max_error(int n, double a,double *d,int problem,vec computed_eigval) {
  double max_error = -1000.0;
  double temp;
  double *analytical_eigval = new double [n];
  computed_eigval = sort(computed_eigval);

  // Beregner analytiske egenvedier FLAGG foreløpig bare buckling beam og quantom dot
  if (problem==0){
    for(int i = 0; i < n; ++i){
      analytical_eigval[i] = analytical_buck(i, n, a, d[i]);
    }
  }
  else if (problem == 1){
    for(int i = 0; i < n; ++i){
      analytical_eigval[i] = analytical_dot(i);
    }
  }
  //beregner maksimal relativ feil
  for(int i = 0; i < n; ++i) {
    temp = fabs((analytical_eigval[i]-computed_eigval[i])/analytical_eigval[i]);
    if (temp > max_error){
      max_error = temp;
    }
  }
  delete[] analytical_eigval;
  return max_error;
}
