#ifndef MAINMPI_H
#define MAINMPI_H

#include <iostream>
#include <mpi.h>
#include <armadillo>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <random>
#include <string>
#include <cmath>

#include "ising.h"
#include "metropolis.h"


using namespace std;
using namespace arma;

void writeExpVals(ofstream&, int&, double&, double&, double&, vec&);

inline bool acceptanceRule(double parameter, double randomNumber){
  return (randomNumber <= parameter);
}

#endif /* MAINMPI.H */
