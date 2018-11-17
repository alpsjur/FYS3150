#include "mainMPI.h"


int main(int argc, char *argv[]){
  // initialising parallellisation
  int myRank, numProcs;
  MPI_Init (&argc, &argv);
  MPI_Comm_size (MPI_COMM_WORLD, &numProcs);
  MPI_Comm_rank (MPI_COMM_WORLD, &myRank);

  // initialising files

  // initialising Ising model
  int gridDimension = 80;
  double couplingParameter = 1;
  double mcCycles = 1e6;
  double equilMC = 25e3;  // found experimentally from studying relaxation time
  double initTemp = 2.0;
  double finalTemp = 2.3;
  double dT = 0.01;
  bool ordered = false;

  IsingModel spinLattice(gridDimension, couplingParameter, ordered);
  spinLattice.initSystem();

  // initialising and opening files
  char *outfilename;
  ofstream outfile;
  if (myRank == 0) {
    outfilename = argv[1];
    outfile.open(outfilename);
  }

  // distributing Monte Carlo cycles
  double totalEquilMC = numProcs*equilMC;  // the equillibrium MC after collecting
  double my_mcCycles = mcCycles/numProcs;
  if ((myRank == numProcs-1) && ((int) mcCycles%numProcs != 0)){
    my_mcCycles += (int) mcCycles%numProcs;
  }

  // parallelising system parameters
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
  // initialising array to hold probability distribution function
  double pdf[2000];
  // running and timing a parallellised Metropolis algorithm
  double timeStartMPI = MPI_Wtime();
  for (double temp = initTemp; temp <= finalTemp; temp += dT){
    expectationValues.zeros(numberOfExpvals);
    collectedExpVals.zeros(numberOfExpvals);
    metropolis(spinLattice, temp, acceptanceRule, my_mcCycles, equilMC, expectationValues, pdf, nodeSeed);
    for(int i = 0; i < numberOfExpvals; ++i){
      MPI_Reduce(&expectationValues[i], &collectedExpVals[i], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    }
    if(myRank == 0){
      writeExpVals(outfile, gridDimension, mcCycles, totalEquilMC, temp, collectedExpVals);
      cout << "temperature = " << temp << " complete" << endl;
    }
  }

  double timeEndMPI = MPI_Wtime();
  if (myRank == 0){
    outfile.close();
    cout << "complete after " << (timeEndMPI - timeStartMPI) << " on " << numProcs << " nodes" << endl;
  }

  MPI_Finalize();
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
