#ifndef ISING_H
#define	ISING_H

#include <armadillo>

using namespace arma;

inline int periodic(int i, int gridDimension){
  return (i+gridDimension) % gridDimension;
}


class IsingModel{
private:
  imat m_spinMatrix;
  double m_couplingParameter;
  int m_rowLength, m_columnLength;
  double m_deltaEnergy, m_deltaMagnetisation;

public:
  double energy, magnetisation;
  IsingModel(imat &spinMatrix, double couplingParameter, double &energy, double &magnetisation)
  : m_spinMatrix(spinMatrix), m_couplingParameter(couplingParameter),
    energy(energy), magnetisation(magnetisation), m_rowLength(spinMatrix.n_rows),
    m_columnLength(spinMatrix.n_cols){}

  void initSystem();
  void flipSpin(int&, int&);
  double calculateDeltaEnergy(int&, int&);
  double calculateDeltaMagnetisation(int&, int&);
  int nRows(){return m_rowLength;}
  int nCols(){return m_columnLength;}
};


#endif /* ISING.H */
