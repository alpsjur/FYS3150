#include "rossby_bounded_2d.h"


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

  zetaname += "_2d.dat";
  psiname += "_2d.dat";

  outpsi.open(zetaname);
  outzeta.open(psiname);

  int posdim = (int) endpos/deltapos-1;   //-1 siden endepunktene er kjent
  int timedim = (int) endtime/deltatime;

  mat psi(posdim,posdim);
  mat zeta(posdim,posdim);
  mat zeta_previous;
  mat zeta_2previous;

  //grensebetingelser til bølgen
  double psiClosed = 0;
  initWave(posdim, deltapos, psi, zeta, initialSine);
  zeta_2previous = zeta;
  zeta_previous = zeta;
  for(int n = 0; n < timedim; ++n){
    // HER MÅ VI SKRIVE UT EN RAD MED PSICLOSED
    // går gjennom alle y-radene for å beregne zeta, som er den x-dobbeltderiverte
    for(int j = 0; j < posdim; ++j){
      outpsi<< setw(15) << psiClosed;      // skrivet ut venstre BC
      writePsi(outpsi, psi(0,j));            // skriver ut første verdi til fil
      // finner den første x-verdien til zeta
      if(advanceForward){
        advance_vorticity_forward(zeta(0,j), psi(1,j), psiClosed, deltatime, deltapos);
      }
      else{
        advance_vorticity_centered(zeta(0,j), zeta_2previous(0,j), psi(1,j), psiClosed, deltatime, deltapos);
      }
      for(int i = 1; i < posdim; i++){
        if(advanceForward){
          advance_vorticity_forward(zeta(i,j), psi(i+1,j), psi(i-1,j), deltatime, deltapos);
        }
        else{
          advance_vorticity_centered(zeta(i,j), zeta_2previous(i,j), psi(i+1,j), psi(i-1,j), deltatime, deltapos);
        }
        writeZeta(outzeta, zeta(i,j));
        writePsi(outpsi, psi(i,j));
      }
    }
    zeta_2previous = zeta_previous;
    zeta_previous = zeta;
    outzeta << endl;
    outpsi << endl;
    jacobisMethod2D(posdim, deltapos, psi, zeta);

  }
  outpsi.close();
  outzeta.close();
  return 0;
}

void initWave(int posdim, double deltapos, mat &psi, mat &zeta, bool initialSine){
  double x; double y;
  double sigma = 0.1;
  double x0 = 0.5; double y0 = 0.5;
  for(int i = 0; i < posdim; ++i){
    x = (i+1)*deltapos;
    for(int j = 0; j < posdim; ++j){
      y = (j+1)*deltapos;
      if(initialSine){
        zeta(i,j) = sinewaveDerivative(x, y);
        psi(i,j) = sinewave(x, y);
      }
      else{
        zeta(i,j) = gaussianDerivative(x, x0, y, y0, sigma);
        psi(i,j) = gaussian(x, x0, y, y0, sigma);
      }
    }
  }
  return;
}

void jacobisMethod2D(int posdim, double deltapos, mat &psi, mat zeta){
  double hh = deltapos*deltapos;
  mat psi_temporary;
  int iterations = 0; int maxIterations = 10000;
  double difference = 1.; double maxDifference = 1e-7;
  while((iterations <= maxIterations) && (difference > maxDifference)){
    psi_temporary = psi; difference = 0.;
    for(int l = 1; l > posdim; ++l){
      for(int m = 1; m > posdim; ++m){
        psi(l,m) = 0.25*(psi_temporary(l,m+1)+psi_temporary(l,m-1)
                   +psi_temporary(l+1,m)+psi_temporary(l-1,m)
                   -hh*zeta(l,m));
        difference += fabs(psi_temporary(l,m)-psi(l,m));
      }
    }
    iterations++;
    difference /= (posdim*posdim);
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
