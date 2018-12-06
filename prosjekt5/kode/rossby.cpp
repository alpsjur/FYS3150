#include "rossby.h"


int main(int argc, char *argv[]){

  ofstream outpsi, outzeta;
  outpsi.open("../data/psi.dat");
  outzeta.open("../data/zeta.dat");

  int timedim = atoi(argv[1]);
  int posdim = atoi(argv[2]);

  double endtime = 20.0;
  double endpos = 1;

  double deltatime = endtime/timedim;
  double deltapos = endpos/posdim;

  // grensebetingelser for vortisiteten og straumfunksjonen
  zeta[0] = 0.0;
  zeta[posdim] = 0.0;


  // initialising a vector of vectors (matrix) to hold a wave for every timestep
  vector<double> psi(posdim);
  vector<double> zeta(posdim);

  // matriseelementer
  /*
  løser likingssett på formen

    [1 0 0 0 ..... 0]   [  v0]   [ 0]
    [a b c 0 ..... 0]   [  v1]   [d0]
    [0 a b c ..... 0]   [  v2]   [d1]
    [0 0 a b ..... 0] * [  v3] = [d2]
    [........... b c]   [....]   [..]
    [0 ......... a b]   [vn+1]   [dn]
    [0 ......... 0 1]   [vn+2]   [ 0]

  der a = c = 1, b = -2

  */
  vector<double> a(posdim+2);
  vector<double> b(posdim+2);
  vector<double> c(posdim+2);
  vector<double> d(posdim+2);
  vector<double> x(posdim);

  initPsi(posdim, psi);
  for(int n = 0; n < timedim; ++n){
    for(int j = 1; j < posdim; ++j){
      //jobbe videre her
      zeta[j] -= (deltatime/deltapos)*(psi[j+1] - psi[j-1]);
      writeZeta(outzeta, zeta[j]);
      writePsi(outpsi, psi[j]);
    }
    outzeta << endl;
    outpsi << endl;
    initMatrixElements(posdim, zeta, x, a, b, c, d);
    forward_sub(posdim, a, b, c, d);
    backward_sub(posdim, psi, b, c, d);
  }
  outpsi.close();
  outzeta.close();
  return 0;
}


//initsaliserer strømfunksjonen
void initPsi(int posdim, vector<double> &psi){
  double x;
  for(int i = 0; i < posdim; ++i){
    x = i*(1./(posdim-1));
    psi[i] = sinewave(x);
    //psi[i] = gaussian(x);
  }
  return;
}


// initialiserer matriseelementene og funksjonsverdiene
void initMatrixElements(int posdim, vector<double> zeta, vector<double> &x, vector<double> &a,
                        vector<double> &b, vector<double> &c, vector<double> &d) {
  const double h = 1.0/(posdim + 1.0);
  const double hh = h*h;
  a[0] = -1.0;
  b[0] = 1.0;
  c[0] = 0.0;
  d[0] = 0.0;
  for(int i = 0; i < posdim; ++i) {
    x[i] = (i + 1)*h;
    a[i+1] = -1.0;
    b[i+1] = 2.0;
    c[i+1] = -1.0;
    d[i+1] = hh*zeta[i];
  }
  a[posdim] = 0;
  b[posdim+1] = 1;
  d[posdim+1] = 0;
  return;
}

// utfører gaussisk eliminasjon for å redusere likningssettet
void forward_sub(int posdim, vector<double> &a, vector<double> &b,
                 vector<double> &c, vector<double> &d) {
  for(int i = 1; i < posdim+2; ++i) {
    const double k = a[i-1]/b[i-1];
    b[i] -= k*c[i-1];
    d[i] -= k*d[i-1];
  }
  return;
}


// løser likningssettet for v-vektoren
void backward_sub(int posdim, vector<double> &psi, vector<double> &b,
                  vector<double> &c, vector<double> &d) {
  for(int i = posdim; i >= 1; --i) {
    psi[i-1] = (d[i] - c[i]*psi[i])/b[i];
  }
  return;
}

void advance_vorticity_forward(double &zeta_forward, double &zeta_center,
                              double &psi_forward, double &psi_backward,
                              double &deltatime, double &deltapos){
  zeta_forward = zeta_center - (deltatime/(2.0*deltapos))*(psi_forward - psi_backward);
  return;
}

void advance_vorticity_centered(double &zeta_forward, double &zeta_backward,
                                double &psi_forward, double &psi_backward,
                                double &deltatime, double &deltapos){
  zeta_forward = zeta_backward - (deltatime/(2.0*deltapos))*(psi_forward - psi_backward);
  return;
}

void writePsi(ofstream &outpsi, double &psivalue){
  outpsi << setiosflags(ios::showpoint | ios::uppercase);
  outpsi << setw(15) << setprecision(8) << psivalue;
  return;
}

void writeZeta(ofstream &outzeta, double &zetavalue){
  outzeta << setiosflags(ios::showpoint | ios::uppercase);
  outzeta << setw(15) << setprecision(8) << zetavalue;
  return;
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
