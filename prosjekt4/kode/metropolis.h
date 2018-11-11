#ifndef METROPOLIS_H
#define METROPOLIS_H

#include <functional>
#include <cmath>
#include "ising.h"

using namespace std;

void metropolis(IsingModel&, double&, function<bool(double, double)>, double&, vec&, long&);

#endif /* METROPOLIS.H */
