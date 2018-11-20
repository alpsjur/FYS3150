#include "main.h"


int main(int argc, char *argv[]){

  // initialising and opening files
  ofstream outfile;
  outfile.open(argv[1]);
  ofstream outPDF;
  //outPDF.open("../data/isingPDF_GRID20_1E6_T24.dat");
  // initialising Ising model
  int gridDimension = 2;
  double couplingParameter = 1;

  double initMC = 25e3;
  double finalMC = 1e6;
  double dMC = 1e4;
  double equilMC = 20e3;  // found experimentally from studying relaxation time

  double initTemp = 1.0;
  double finalTemp = 1.0;
  double dT = 0.01;
  bool ordered = false;

  IsingModel spinLattice(gridDimension, couplingParameter, ordered);
  spinLattice.initSystem();

  //initialising random seed
  long nodeSeed = -1;

  // initialising vectors to hold expectation values obtained from MC
  int numberOfExpvals = 5;
  vec expectationValues(numberOfExpvals, fill::zeros);
  // initialising array to hold probability distribution function
  double pdf[2000] = {0};
  // running the metropolis algorithm for different mcs and temperatures
  for(double mcCycles = initMC; mcCycles <= finalMC; mcCycles += dMC){
    for (double temp = initTemp; temp <= finalTemp; temp += dT){
      expectationValues.zeros(numberOfExpvals);
      metropolis(spinLattice, temp, acceptanceRule, mcCycles, equilMC, expectationValues, pdf, nodeSeed);
      //writeExpVals(outfile, gridDimension, mcCycles, equilMC, temp, expectationValues);
      //writePDF(outPDF, pdf, gridDimension);
      writeRelativeError(outfile, mcCycles, equilMC, expectationValues);
    }
    //writeMCexpvals(outfile, mcCycles, gridDimension, expectationValues);
  }
  outfile.close();
  outPDF.close();

  return 0;
}


void writeExpVals(ofstream &outfile, int &gridDimension, double &mcCycles, double &equilMC, double &temperature, vec &expectationValues){
  double gridSize = double(gridDimension*gridDimension);
  double norm = 1.0/((double) (mcCycles - equilMC));  // divided by  number of cycles
  double energy_expvals = expectationValues(0)*norm;
  double energy2_expvals = expectationValues(1)*norm;
  double magnetisation_expvals = expectationValues(2)*norm;
  double magnetisation2_expvals = expectationValues(3)*norm;
  double magnetisationAbs_expval = expectationValues(4)*norm;
  // all expectation values are per spin, divide by 1/NSpins/NSpins
  double energyVariance = (energy2_expvals - energy_expvals*energy_expvals)/gridSize;
  double magnetisationVariance = (magnetisation2_expvals - magnetisationAbs_expval*magnetisationAbs_expval)/gridSize;
  outfile << setiosflags(ios::showpoint | ios::uppercase);
  outfile << setw(15) << setprecision(8) << temperature;
  outfile << setw(15) << setprecision(8) << energy_expvals/gridSize;
  outfile << setw(15) << setprecision(8) << energyVariance/(temperature*temperature);
  outfile << setw(15) << setprecision(8) << magnetisation_expvals/gridSize;
  outfile << setw(15) << setprecision(8) << magnetisationVariance/temperature;
  outfile << setw(15) << setprecision(8) << magnetisationAbs_expval/gridSize << endl;

  return;
}


void writePDF(ofstream &outfile, double *pdf, int &gridDimension){
  double gridSize = double(gridDimension*gridDimension);
  double energy = 0;
  for(int i = 0; i < 2000; ++i){
    if(pdf[i] != 0){
      if(i > 1000){
        energy = (double)  -(i - 1000.0)/gridSize;
      }
      else {
        energy = (double) i/gridSize;
      }
      outfile << setw(15) << setprecision(8) << (double) energy;
      outfile << setw(15) << setprecision(8) << pdf[i] << endl;
    }
  }
  return;
}


void writeMCexpvals(ofstream &outfile, double &mcCycles, int &gridDimension, vec &expectationValues){
  double gridSize = double(gridDimension*gridDimension);
  double norm = 1.0/(double) mcCycles;  // divided by  number of cycles
  double energy_expvals = expectationValues(0)*norm;
  double magnetisationAbs_expval = expectationValues(4)*norm;
  outfile << setiosflags(ios::showpoint | ios::uppercase);
  outfile << setw(15) << setprecision(8) << mcCycles;
  outfile << setw(15) << setprecision(8) << energy_expvals/gridSize;
  outfile << setw(15) << setprecision(8) << magnetisationAbs_expval/gridSize << endl;

  return;
}

void writeRelativeError(ofstream &outfile, double &mcCycles, double& equilMC, vec &expectationValues){
  double couplingParameter = 1.0;
  double temperature = 1.0;
  double boltzmannsConstant = 1.0;
  double beta = 1.0/(temperature*boltzmannsConstant);
  double betaJ = 8.0*beta*couplingParameter;
  double meanEnergy = - 8.0*couplingParameter*sinh(betaJ)
                      /(cosh(betaJ) + 3);
  double meanMagnetisation = (2.0*exp(-betaJ) + 4);

  double norm = 1.0/((double) (mcCycles - equilMC));
  double energy_expvals = expectationValues(0)*norm;
  double magnetisationAbs_expval = expectationValues(4)*norm;

  outfile << setiosflags(ios::showpoint | ios::uppercase);
  outfile << setw(15) << setprecision(8) << mcCycles;
  outfile << setw(15) << setprecision(8) << fabs(meanEnergy - energy_expvals)/fabs(meanEnergy);
  outfile << setw(15) << setprecision(8) << fabs(meanMagnetisation - magnetisationAbs_expval)/fabs(meanMagnetisation) << endl;

  return;
}
