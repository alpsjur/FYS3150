#include "metropolis.h"


void metropolis(IsingModel &spinLattice, double &temperature, function<bool(double,
                double)> acceptanceRule, double &mcCycles, vec &expectationValues,
                long &nodeSeed)
  {
  // calling the Mersienne algo
  mt19937_64 gen(nodeSeed);
  // Set up the uniform distribution for x \in [[0, 1]
  uniform_real_distribution<double> RNG(0.0,1.0);

  // initialising possible energy differences
  vec energyDifferences = zeros<mat>(17);
  for( int dE =-8; dE <= 8; dE += 4) {energyDifferences(dE+8) = exp(-dE/temperature);}

  int rowRandom, columnRandom, deltaEnergy, deltaMagnetisation;
  for(int cycles = 0; cycles < mcCycles; ++cycles){
    for(int row = 0; row < spinLattice.nRows(); ++row){
      for(int column = 0; column < spinLattice.nCols(); ++column){
        rowRandom = (int) (RNG(gen)*(double)spinLattice.nRows());
        columnRandom = (int) (RNG(gen)*(double)spinLattice.nCols());
        deltaEnergy = spinLattice.calculateDeltaEnergy(rowRandom, columnRandom);

        if(acceptanceRule(energyDifferences(deltaEnergy + 8), RNG(gen))){
          spinLattice.flipSpin(rowRandom, columnRandom);
          deltaMagnetisation = spinLattice.calculateDeltaMagnetisation(rowRandom, columnRandom);
          spinLattice.magnetisation += deltaMagnetisation;
          spinLattice.energy += deltaEnergy;
        }
      }
    }
    expectationValues(0) += spinLattice.energy;
    expectationValues(1) += spinLattice.energy * spinLattice.energy;
    expectationValues(2) += spinLattice.magnetisation;
    expectationValues(3) += spinLattice.magnetisation * spinLattice.magnetisation;
    expectationValues(4) += fabs(spinLattice.magnetisation);
  }
  return;
}
