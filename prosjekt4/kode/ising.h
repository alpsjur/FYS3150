#ifndef ISING_H
#define	ISING_H

#include <armadillo>

using namespace arma;

inline int periodic(int i, int gridDimension){
  return (i+gridDimension) % gridDimension;
}


class IsingModel{
  // class that initialises an Ising lattice, as well as provide functions for simulation
private:
  imat m_spinMatrix;
  double m_couplingParameter;
  int m_rowLength, m_columnLength;
  double m_deltaEnergy, m_deltaMagnetisation;
  bool m_ordered;

public:
  double energy, magnetisation;
  IsingModel(int &latticeDimension, double &couplingParameter, bool &ordered)
  : m_couplingParameter(couplingParameter), m_ordered(ordered){
      m_spinMatrix.zeros(latticeDimension, latticeDimension);
      m_rowLength = m_spinMatrix.n_rows;
      m_columnLength = m_spinMatrix.n_cols;
    }

  void initSystem();
  void flipSpin(int&, int&);
  double calculateDeltaEnergy(int&, int&);
  double calculateDeltaMagnetisation(int&, int&);
  int nRows(){return m_rowLength;}
  int nCols(){return m_columnLength;}
};


#endif /* ISING.H */
