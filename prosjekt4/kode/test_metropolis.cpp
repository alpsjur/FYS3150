#include <mpi.h>
#include <cmath>
#include <armadillo>
#include <vector>
#include <functional>

#include "catch.hpp"
#include "ising.h"
#include "metropolis.h"

using namespace std;
using namespace arma;

inline bool acceptanceRule(double parameter, double randomNumber){
  return (randomNumber <= parameter);
}


TEST_CASE("ANALYTICAL 2X2 ISING LATTICE"){
  // analytical expressions
  double couplingParameter = 1;
  double temperature = 1;
  double boltzmannsConstant = 1;
  double beta = temperature*boltzmannsConstant;
  double betaJ = 8*beta*couplingParameter;
  double meanEnergy = - 4.0*couplingParameter*sinh(betaJ)
                      /(cosh(betaJ) + 6);
  double heatCapacitance = 32.0*couplingParameter*couplingParameter/(beta*temperature)
                          * (sinh(betaJ)*sinh(betaJ) - cosh(betaJ)*cosh(betaJ) - 3.0*cosh(betaJ))
                          / (cosh(betaJ) + 3)*(cosh(betaJ) + 3);
  double meanMagnetisation = 0;


  double energy, magnetisation;
  int gridDimension = 2;

  imat spinMatrix(gridDimension, gridDimension, fill::ones);
  IsingModel spinLattice(spinMatrix, couplingParameter, energy, magnetisation);
  spinLattice.initSystem();


  double mcCycles = 1e6;
  vec expectationValues(5, fill::zeros);

  metropolis(spinLattice, temperature, acceptanceRule, mcCycles, expectationValues);

  double norm = 1.0/((double) (mcCycles));
  double gridSize = double(gridDimension*gridDimension);
  double energy_expvals = expectationValues(0)*norm/gridSize;
  double magnetisation_expvals = expectationValues(2)*norm/gridSize;

  REQUIRE(energy_expvals == meanEnergy);
  REQUIRE(magnetisation_expvals == meanMagnetisation);
}
