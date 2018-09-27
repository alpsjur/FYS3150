#include "catch.hpp"
#include "jacobis_method.hpp"

TEST_CASE("Buckling beam : Sjekker om vi finner største verdi i matrisen"){
    int n = 3;
    int problem = 0;

    double rhomin = 0.0;
    double rhomax = 1.0;
    double a, omega_r;

    int k, l;

    double *d = new double [n];
    double *analytical_eigval = new double [n];

    mat A(n,n, fill::zeros);

    //initialize matrices and vector

    initialize(A, d, a, rhomin, rhomax, omega_r, problem, n);
    delete[] analytical_eigval, d;
    //find maximum matrix element
    double largest = find_largest(A, k, l, n);

    REQUIRE(k==0);
    REQUIRE(l==1);
    REQUIRE(largest==Approx(9));
}


TEST_CASE("Buckling beam : analytical eigval == numerical eigval"){
  int n = 3;
  int problem = 0;

  double rhomin = 0.0;
  double rhomax = 1.0;

  double a, omega_r;

  mat A(n, n, fill::zeros);
  mat S(n, n, fill::eye);
  vec jacobi_eigval(n);

  double *d = new double [n];
  double *analytical_eigval = new double [n];


  initialize(A, d, a, rhomin, rhomax, omega_r, problem, n);

  for(int i = 0; i < n; ++i){
    analytical_eigval[i] = analytical_buck(i, n, a, d[i]);
  }
  delete[] d;

  int iterations;
  jacobi(n, iterations, A, S, jacobi_eigval);
  jacobi_eigval = sort(jacobi_eigval);

  for(int i = 0; i<n; ++i){
    REQUIRE(analytical_eigval[i] == Approx(jacobi_eigval[i]));
  }
  delete[] analytical_eigval;
}


TEST_CASE("Quantum dot : Sjekker om vi finner største verdi i matrisen"){
    int n = 3;
    int problem = 1;

    double rhomin = 0.0;
    double rhomax = 10.0;
    double a, omega_r;

    int k, l;

    double *d = new double [n];
    double *analytical_eigval = new double [n];

    mat A(n,n, fill::zeros);

    //initialize matrices and vector

    initialize(A, d, a, rhomin, rhomax, omega_r, problem, n);
    delete[] analytical_eigval, d;
    //find maximum matrix element
    double largest = find_largest(A, k, l, n);

    REQUIRE(k==0);
    REQUIRE(l==1);
    REQUIRE(largest==Approx(0.09));
}

TEST_CASE("Quantun dot: analytiske eigval == numeriske eigval"){
  int n = 180;
  int problem = 1;

  double rhomin = 0.0;
  double rhomax = 10;

  double a, omega_r;

  mat A(n, n, fill::zeros);
  mat S(n, n, fill::eye);
  vec jacobi_eigval(n);

  double *d = new double [n];
  double *analytical_eigval = new double [n];


  initialize(A, d, a, rhomin, rhomax, omega_r, problem, n);

  analytical_eigval[0] = 3.0;
  for(int i = 1; i < n; ++i){
    analytical_eigval[i] = analytical_eigval[i-1] +  4.0;
  }
  delete[] d;

  int iterations;
  jacobi(n, iterations, A, S, jacobi_eigval);
  jacobi_eigval = sort(jacobi_eigval);
  for(int i = 0; i<3; ++i){
    REQUIRE(analytical_eigval[i] == Approx(jacobi_eigval[i]).epsilon(0.001));
  }
  delete[] analytical_eigval;
}
