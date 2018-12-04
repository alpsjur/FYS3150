#include "rossby.h"

int main(int argc, char *argv){

  ofstream outfile;
  outfile.open(argv[3]);

  int timedim = atoi(argv[1]);
  int posdim = atoi(argv[2]);


  // initialising a vector of vectors (matrix) to hold a wave for every timestep
  vector<double> wave;

  for(int n = 0; n < timedim; ++n){
    for(int j = 0; j < posdim; ++j){

    }
  }
  outfile.close()
  return 0;
}


// initialiserer matriseelementene og funksjonsverdiene
void initMatrixElements(int posdim, double *x, double *a, double *b, double *c, double *d) {
  const double h = 1.0/(n + 1.0);
  const double hh = h*h;

  for(int i = 0; i < dim; ++i) {
    x[i] = (i + 1)*h;
    a[i] = -1.0;
    b[i] = 2.0;
    c[i] = -1.0;
    d[i] = hh*f(x[i]);
  }
}

// utfører gaussisk eliminasjon for å redusere likningssettet
void forward_sub(int dim, double *a, double *b, double *c, double *d) {
  for(int i = 1; i < dim; ++i) {
    const double k = a[i-1]/b[i-1];
    b[i] -= k*c[i-1];
    d[i] -= k*d[i-1];
  }
}


// løser likningssettet for v-vektoren
void backward_sub(int dim, double *v, double *b, double *c, double *d) {
  for(int i = dim-2; i >= 0; --i) {
    v[i] = (d[i] - c[i]*v[i+1])/b[i];
  }
}

void advance_vorticity_forward(double &xi_forward, double &xi_center,
                              double &psi_forward, double &psi_backward,
                              double &deltatime, double &deltapos){
  xi_forward = xi_center - (deltatime/(2.0*deltapos))*(psi_forward - psi_backward);
  return;
}

void advance_vorticity_centered(double &xi_forward, double &xi_backward,
                                double &psi_forward, double &psi_backward,
                                double &deltatime, double &deltapos){
  xi_forward = xi_backward - (deltatime/(2.0*deltapos))*(psi_forward - psi_backward);
  return;
}

void calculate_error(int n, double *v, double *x, double *eps) {
  for(int i = 0; i < n; ++i) {
    eps[i] = log10(fabs((exact(x[i]) - v[i])/exact(x[i])));
  }
}

void centered_difference(double &derivative, double &forward, double &backward, double &step){
  derivative = (forward - backward)/(2.0*step);
  return;
}

void centered_difference2(double &derivative, double &forward, double &center, double &backward, double &step){
  derivative = (forward - 2.0*center + backward)/(step*step);
  return;
}

void forward_difference(double &derivative, double &forward, double &center, double &step){
  derivative = (forward - center)/step;
  return;
}
