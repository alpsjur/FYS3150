#ifndef MAIN_H
#define MAIN_H

#include <iostream>
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
void writePDF(ofstream&, double*, int&, double&, double&, vec&);
void writeMCexpvals(ofstream&, double&, int&, vec&);
void writeRelativeError(ofstream&, double&, double&, vec&);
void writeAcceptanceRate(ofstream&, double&, double&, int&);

inline bool acceptanceRule(double parameter, double randomNumber){
  return (randomNumber <= parameter);
}

#endif /* MAIN.H */
