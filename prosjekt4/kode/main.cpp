#include "main.h"


int main(int argc, char *argv[]){
  // initialising parallellisation
  int myRank, numProcs;
  MPI_Init (&argc, &argv);
  MPI_Comm_size (MPI_COMM_WORLD, &numProcs);
  MPI_Comm_rank (MPI_COMM_WORLD, &myRank);

  // initialising files

  // initialising Ising model
  double energy, magnetisation;
  int gridDimension = 2;
  double couplingParameter = 1;

  imat spinMatrix(gridDimension, gridDimension, fill::ones);
  IsingModel spinLattice(spinMatrix, couplingParameter, energy, magnetisation);
  spinLattice.initSystem();

  // initialising and opening files
  char *outfilename;
  ofstream outfile;
  if (myRank == 0) {
    outfilename = argv[1];
    outfile.open(outfilename);
  }

  // initialising and distributing Monte Carlo cycles
  double mcCycles = 10000;
  double my_mcCycles = mcCycles/numProcs;
  if ((myRank == numProcs-1) && ((int) mcCycles%numProcs != 0)){
    my_mcCycles += (int) mcCycles%numProcs;
  }

  // initialising and parallelisising temperature parameters
  double initTemp = 1.0;
  double finalTemp = 2.5;
  double dT = 0.5;
  MPI_Bcast (&gridDimension, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast (&initTemp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast (&finalTemp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast (&dT, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  // distributing different random seeds for each node
  long nodeSeed = -1 - myRank;

  // initialising vectors to hold expectation values obtained from MC
  int numberOfExpvals = 5;
  vec expectationValues(numberOfExpvals, fill::zeros);
  vec collectedExpVals(numberOfExpvals, fill::zeros);

  // running and timing a parallellised Metropolis algorithm
  double timeStartMPI = MPI_Wtime();
  for (double temp = initTemp; temp <= finalTemp; temp += dT){
    vec expectationValues(5, fill::zeros);
    metropolis(spinLattice, temp, acceptanceRule, mcCycles, expectationValues, nodeSeed);
    for(int i = 0; i < numberOfExpvals; ++i){
      MPI_Reduce(&expectationValues[i], &collectedExpVals[i], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    }
    if(myRank == 0){
      writeToFile(outfile, gridDimension, mcCycles, temp, collectedExpVals);
    }
  }

  double timeEndMPI = MPI_Wtime();
  if (myRank == 0){
    outfile.close();
    cout << "complete after " << (timeEndMPI - timeStartMPI) << " on " << numProcs << " nodes";
  }

  MPI_Finalize();
  return 0;
}


void writeToFile(ofstream &outfile, int gridDimension, int mcCycles, double temperature, vec& expectationValues){
  double gridSize = double(gridDimension*gridDimension);
  double norm = 1.0/((double) (mcCycles));  // divided by  number of cycles
  double energy_expvals = expectationValues(0)*norm;
  double energy2_expvals = expectationValues(1)*norm;
  double magnetisation_expvals = expectationValues(2)*norm;
  double magnetisation2_expvals = expectationValues(3)*norm;
  double magnetisationAbs_expval = expectationValues(4)*norm;
  // all expectation values are per spin, divide by 1/NSpins/NSpins
  double energyVariance = (energy_expvals - energy_expvals*energy_expvals)/gridSize;
  double magnetisationVariance = (magnetisation2_expvals - magnetisationAbs_expval*magnetisationAbs_expval)/gridSize;
  outfile << setiosflags(ios::showpoint | ios::uppercase);
  outfile << setw(15) << setprecision(8) << temperature;
  outfile << setw(15) << setprecision(8) << energy_expvals/gridSize;
  outfile << setw(15) << setprecision(8) << energyVariance/temperature/temperature;
  outfile << setw(15) << setprecision(8) << magnetisation_expvals/gridSize;
  outfile << setw(15) << setprecision(8) << magnetisationVariance/temperature;
  outfile << setw(15) << setprecision(8) << magnetisationAbs_expval/gridSize << endl;
}
