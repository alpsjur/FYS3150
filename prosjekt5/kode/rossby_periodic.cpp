#include "rossby_periodic.h"


// reading n power from command line
int main(int argc, char *argv[]) {
  ofstream outpsi, outzeta;
  string zetaname = "../data/zeta_periodic";
  string psiname = "../data/psi_periodic";

  double deltapos = atof(argv[1]);
  double deltatime = atof(argv[2]);
  double endtime = atof(argv[3]);
  double endpos = 1.0;

  bool initialSine;
  if(atof(argv[4])==0){
    initialSine = true;
    zetaname += "_sine";
    psiname += "_sine";
  }
  else{
    initialSine = false;
    zetaname += "_gaussian";
    psiname += "_gaussian";
  }

  bool advanceForward;
  if(atof(argv[5])==0){
    advanceForward = true;
    zetaname += "_forward";
    psiname += "_forward";
  }
  else{
    advanceForward = false;
    zetaname += "_centered";
    psiname += "_centered";
  }

  zetaname += ".dat";
  psiname += ".dat";

  outpsi.open(psiname);
  outzeta.open(zetaname);

  int posdim = (int) (endpos/deltapos);
  int timedim = (int) endtime/deltatime;

  vector<double> psi(posdim);
  vector<double> zeta(posdim);
  vector<double> zeta_previous;
  vector<double> zeta_2previous;

  initWave(posdim, deltapos, psi, zeta, initialSine);
  zeta_2previous = zeta;
  zeta_previous = zeta;
  for(int n = 0; n < timedim; ++n){
    // finner den første x-verdien til zeta
    for(int j = 0; j < posdim; ++j){
      if(advanceForward){
        advance_vorticity_forward(zeta[j], psi[periodic(j+1, posdim)], psi[periodic(j-1, posdim)], deltatime, deltapos);
      }
      else{
        advance_vorticity_centered(zeta[j], zeta_2previous[j], psi[periodic(j+1, posdim)], psi[periodic(j-1, posdim)], deltatime, deltapos);
      }
      writeZeta(outzeta, zeta[j]);
      writePsi(outpsi, psi[j]);
    }
    zeta_2previous = zeta_previous;
    zeta_previous = zeta;
    outpsi << setw(15) << psi[0];    // skriver ut høyre BC
    outzeta << setw(15) << zeta[0];    // skriver ut høyre BC
    outzeta << endl;
    outpsi << endl;
    jacobisMethod(posdim, deltapos, psi, zeta);
  }
  outpsi.close();
  outzeta.close();
  return 0;
}


void initWave(int posdim, double deltapos, vector<double> &psi, vector<double> &zeta, bool initialSine){
  double x;
  //double h = 1.0/(posdim + 1.0);
  double sigma = 0.1; double x0 = 0.5;
  for(int j = 0; j < posdim; ++j){
    x = j*deltapos;
    if(initialSine){
      zeta[j] = sinewaveDerivative(x);
      psi[j] = sinewave(x);
    }
    else{
      zeta[j] = gaussianDerivative(x, x0, sigma);
      psi[j] = gaussian(x, x0, sigma);
    }
  }
  return;
}
void jacobisMethod(int posdim, double deltapos, vector<double> &psi, vector<double> zeta){
  double hh = deltapos*deltapos;
  vector<double> psi_temporary;
  int iterations = 0; int maxIterations = 10000;
  double difference = 1.; double maxDifference = 1e-6;
  while((iterations <= maxIterations) && (difference > maxDifference)){
    psi_temporary = psi; difference = 0.;
    for(int j = 0; j < posdim; ++j){
      psi[j] = 0.5*(psi_temporary[periodic(j+1,posdim)]+psi_temporary[periodic(j-1,posdim)]-zeta[j]*hh);
      difference += fabs(psi_temporary[j]-psi[j]);
    }
    iterations++;
    difference /= posdim;
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


void writePsi(ofstream &outpsi, double psivalue){
  outpsi << setiosflags(ios::showpoint | ios::uppercase);
  outpsi << setw(15) << setprecision(8) << psivalue;
  return;
}

void writeZeta(ofstream &outzeta, double zetavalue){
  outzeta << setiosflags(ios::showpoint | ios::uppercase);
  outzeta << setw(15) << setprecision(8) << zetavalue;
  return;
}
