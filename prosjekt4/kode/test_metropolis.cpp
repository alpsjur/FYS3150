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
  double couplingParameter = 1.0;
  double temperature = 1.0;
  double boltzmannsConstant = 1.0;
  double beta = 1.0/(temperature*boltzmannsConstant);
  double betaJ = 8.0*beta*couplingParameter;
  double meanEnergy = - 8.0*couplingParameter*sinh(betaJ)
                      /(cosh(betaJ) + 3);
  double heatCapacity = 64.0*couplingParameter*couplingParameter/(beta*temperature)
                        * (1.0 + 3.0*cosh(betaJ))/((cosh(betaJ) + 3)*(cosh(betaJ) + 3));
  double meanMagnetisation = (2.0*exp(betaJ) + 4)/(3 + cosh(betaJ));
  double susceptibility = (8*(exp(betaJ) + 1)/(3 + cosh(betaJ)) - meanMagnetisation*meanMagnetisation)*beta;


  int gridDimension = 2;
  bool ordered = false;

  IsingModel spinLattice(gridDimension, couplingParameter, ordered);
  spinLattice.initSystem();

  double equilMC = 25e3;
  double mcCycles = 1e6;

  long nodeSeed = -1;

  double pdf[2000];
  int acceptanceCounter;
  // initialising vectors to hold expectation values obtained from MC
  int numberOfExpvals = 5;
  vec expectationValues(numberOfExpvals, fill::zeros);

  // running metropolis algorithm
  metropolis(spinLattice, temperature, acceptanceRule, mcCycles, equilMC, expectationValues, pdf, acceptanceCounter, nodeSeed);

  double norm = 1.0/((double) (mcCycles - equilMC));
  double energy_expval = expectationValues(0)*norm;
  double energy2_expval = expectationValues(1)*norm;
  double mcHeatCapacity = (energy2_expval - energy_expval*energy_expval)*(beta/temperature);
  double magnetisation_expval = expectationValues(4)*norm;
  double magnetisation2_expvals = expectationValues(3)*norm;
  double susceptibility_num = (magnetisation2_expvals - magnetisation_expval*magnetisation_expval)/temperature;

  REQUIRE(energy_expval == Approx(meanEnergy).epsilon(0.01));
  REQUIRE(magnetisation_expval == Approx(meanMagnetisation).epsilon(0.01));
  REQUIRE(mcHeatCapacity == Approx(heatCapacity).epsilon(0.01));
  REQUIRE(susceptibility_num == Approx(susceptibility).epsilon(0.01));
}
