#include "rossby_periodic.h"


// reading n power from command line
int main(int argc, char *argv[]) {
  ofstream outpsi, outzeta;

  double deltapos = atof(argv[1]);
  double deltatime = atof(argv[2]);
  double endtime = atof(argv[3]);
  double endpos = 1.0;

  bool initialSine;
  if(atof(argv[4])==0){
    initialSine = true;
  }
  else{
    initialSine = false;
  }

  bool advanceForward;
  if(atof(argv[5])==0){
    advanceForward = true;
    outpsi.open("../data/psi_periodic_forward.dat");
    outzeta.open("../data/zeta_periodic_forward.dat");
  }
  else{
    advanceForward = false;
    outpsi.open("../data/psi_periodic_centered.dat");
    outzeta.open("../data/zeta_periodic_centered.dat");
  }

  int posdim = (int) endpos/deltapos;
  int timedim = (int) endtime/deltatime;

  vector<double> zeta(posdim);
  vector<double> zeta_previous;
  vector<double> zeta_2previous;

  double psiClosed = 0;

  // declaring matrices
  mat A(posdim, posdim, fill::zeros);    //
  mat L;                       // lower triangular matrix
  mat U;                       // upper triangular matrix

  // declaring vectors
  vec f(posdim);                    // column vector containing function values
  vec y(posdim);                    // solution of Ly = f
  vec psi(posdim);                    // solution of Ux = y


  initWave(posdim, psi, zeta, initialSine);
  if(!advanceForward){
    zeta_2previous = zeta;
    zeta_previous = zeta;
  }
  for(int n = 0; n < timedim; ++n){
    // finner den første x-verdien til zeta
    for(int j = 0; j < posdim; ++j){
      if(advanceForward){
        advance_vorticity_forward(zeta[j], psi[periodic(j+1, posdim)], psi[periodic(j-1, posdim)], deltatime, deltapos);
      }
      else{
        advance_vorticity_centered(zeta[j], zeta_2previous[j], psi[periodic(j+1, posdim)], psi[periodic(j-1, posdim)], deltatime, deltapos);
        zeta_2previous = zeta_previous;
        zeta_previous = zeta;
      }
      writeZeta(outzeta, zeta[j]);
      writePsi(outpsi, psi[j]);
    }
    outzeta << endl;
    outpsi << endl;
    initialise(posdim, A, f, zeta);    // initialising A with tridiagonal values
    lu(L, U, A);                 // performing LU-decomposition on A
    solve(y, L, f);     // solving for y indicating that L is triangular
    solve(psi, U, y);
  }
  outpsi.close();
  outzeta.close();
  return 0;
}


void initWave(int posdim, vec &psi, vector<double> &zeta, bool initialSine){
  double x;
  double h = 1.0/(posdim + 1.0);
  double sigma = 0.1;
  double sigma2 = sigma*sigma;
  for(int j = 0; j < posdim; ++j){
    x = (j + 1)*h;
    if(initialSine){
      zeta[j] = -16.0*pi*pi*sinewave(x);
      psi[j] = sinewave(x);
    }
    else{
      zeta[j] = -2.0*gaussian(x, sigma)*(sigma2 - 2.0*x*x)/(sigma2*sigma2);
      psi[j] = gaussian(x, sigma);
    }
  }
  return;
}

void initialise(int posdim, mat &A, vec &f, vector<double> zeta) {
  const double h = 1.0/(posdim + 1.0);
  const double hh = h*h;
  A(0,0) = 2.0;
  A(0, posdim-1) = -1;
  A(posdim-1, 0) = -1;
  for (int j=1; j < posdim; ++j) {
    f[j] = hh*zeta[j];
    A(j,j) = -2.0;
    A(j-1,j) = 1.0;
    A(j, j-1) = 1.0;
  }
  return;
}

void jacobisMethod2D(int posdim, double deltapos, mat &psi, mat zeta){
  double hh = deltapos*deltapos;
  mat psi_temporary;
  int iterations = 0; int maxIterations = 10000;
  double difference = 1.; double maxDifference = 1e-5;
  while((iterations <= maxIterations) && (difference > maxDifference)){
    psi_temporary = psi; difference = 0.;
    for(int l = 1; l > posdim; ++l){
      for(int m = 1; m > posdim; ++m){
        psi(l,m) = 0.25*(psi_temporary(l,m+1)+psi_temporary(l,m-1)
                   +psi_temporary(l+1,m)+psi_temporary(l-1,m)
                   -hh*zeta(l,m));
      }
    }
    iterations++;
    difference /= pow(posdim,2.0);
  }
  return;
}

void advance_vorticity_forward(double &zeta_forward, double psi_forward,
                              double psi_backward, double deltatime,
                              double deltapos){
  zeta_forward -= (deltatime/(2.0*deltapos))*(psi_forward - psi_backward);
  return;
}

void advance_vorticity_centered(double &zeta_forward, double zeta_backward,
                                double psi_forward, double psi_backward,
                                double deltatime, double deltapos){
  zeta_forward = (deltatime/deltapos)*(psi_backward-psi_forward) + zeta_backward;
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
