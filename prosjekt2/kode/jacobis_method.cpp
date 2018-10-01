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
  jacobi_eigval = A.diag();                //vektor med egenverdiene
  //jacobi_eigval = sort(jacobi_eigval);     //sorterer i stigende rekkefølge
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
    if(problem == 1 or problem == 2){
      d[i] += omega_r*V;
      if(interact == 1 and problem == 2){
        d[i] += 1.0/rho;
      }
    }
  }
  //setter diagonalelementene til a og elementene direkte over og under til d
  A(0,0) = d[0];
  for (int i=1; i < n; ++i) {
    //beregner analytiske egenverdier i stigende rekkefølge
    A(i,i) = d[i];
    A(i-1,i) = a;
    A(i, i-1) = a;
  }
  return;
}

double find_largest(mat A, int &k, int &l, int n) {
  //finner indeksene til det største elementet i matrisen
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
    double sik;
    sik = S(i,k);
    S(i,k) = c*sik - s*S(i,l);
    S(i,l) = c*S(i,l) + s*sik;
  }
  return;
}


void write_log(int n, int iterations, double arma_time, double jacobi_time,
                double max_error){
  ofstream logg;
  logg.open("../data/jacobi_log.dat", fstream::app);
  logg << n << ' ' << iterations << ' ' << jacobi_time << ' ' << arma_time
       << ' ' << max_error << endl;
  logg.close();
}

void write_eig(int n, int prob,vec jacobi_eig, mat S){
  ofstream file;
  file.open("../data/jacobi_eig"+ to_string(prob)+ '_' + to_string(n) + ".dat");
  for (int i=0;i<n;++i){
    file << jacobi_eig(i) << ' ' << S.col(i).t();
  }
  file.close();
}

double calculate_max_error(int n, vec computed_eigval, double *analytical_eigval) {
  double max_error = -1000.0;
  double temp;
  computed_eigval = sort(computed_eigval);
  for(int i = 0; i < n; ++i) {
    temp = fabs((analytical_eigval[i]-computed_eigval[i])/analytical_eigval[i]);
    if (temp > max_error){
      max_error = temp;
    }
  }
  return max_error;
}
