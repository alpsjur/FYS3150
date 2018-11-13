#include "metropolis.h"


void metropolis(IsingModel &spinLattice, double &temperature, function<bool(double,
                double)> acceptanceRule, double &mcCycles, double &equilMC,
                vec &expectationValues, long &nodeSeed)
  {

  // initialise possible energy changes
  vec energyChanges = zeros<mat>(17);
  for(int dE =-8; dE <= 8; dE += 4) {energyChanges(dE+8) = exp(-dE/temperature);}

  int rowRandom, columnRandom, deltaEnergy, deltaMagnetisation;
  for(int cycles = 0; cycles < mcCycles; ++cycles){
    for(int row = 0; row < spinLattice.nRows(); ++row){
      for(int column = 0; column < spinLattice.nCols(); ++column){
        rowRandom = (int) (ran2(&nodeSeed)*(double)spinLattice.nRows());
        columnRandom = (int) (ran2(&nodeSeed)*(double)spinLattice.nCols());
        deltaEnergy = spinLattice.calculateDeltaEnergy(rowRandom, columnRandom);

        if(acceptanceRule(energyChanges(deltaEnergy + 8), ran2(&nodeSeed))){
          spinLattice.flipSpin(rowRandom, columnRandom);
          deltaMagnetisation = spinLattice.calculateDeltaMagnetisation(rowRandom, columnRandom);
          spinLattice.magnetisation += deltaMagnetisation;
          spinLattice.energy += deltaEnergy;
        }
      }
    }
    if(cycles > equilMC){
      expectationValues(0) += spinLattice.energy;
      expectationValues(1) += spinLattice.energy * spinLattice.energy;
      expectationValues(2) += spinLattice.magnetisation;
      expectationValues(3) += spinLattice.magnetisation * spinLattice.magnetisation;
      expectationValues(4) += fabs(spinLattice.magnetisation);
    }
  }
  return;
}


#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

double ran2(long *idum)
{
  int            j;
  long           k;
  static long    idum2 = 123456789;
  static long    iy=0;
  static long    iv[NTAB];
  double         temp;

  if(*idum <= 0) {
    if(-(*idum) < 1) *idum = 1;
    else             *idum = -(*idum);
    idum2 = (*idum);
    for(j = NTAB + 7; j >= 0; j--) {
      k     = (*idum)/IQ1;
      *idum = IA1*(*idum - k*IQ1) - k*IR1;
      if(*idum < 0) *idum +=  IM1;
      if(j < NTAB)  iv[j]  = *idum;
    }
    iy=iv[0];
  }
  k     = (*idum)/IQ1;
  *idum = IA1*(*idum - k*IQ1) - k*IR1;
  if(*idum < 0) *idum += IM1;
  k     = idum2/IQ2;
  idum2 = IA2*(idum2 - k*IQ2) - k*IR2;
  if(idum2 < 0) idum2 += IM2;
  j     = iy/NDIV;
  iy    = iv[j] - idum2;
  iv[j] = *idum;
  if(iy < 1) iy += IMM1;
  if((temp = AM*iy) > RNMX) return RNMX;
  else return temp;
}
#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX
