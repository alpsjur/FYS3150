#include "ising.h"

void IsingModel::initSystem(){
  magnetisation = 0;
  energy = 0;
  for(int row = 0; row < m_rowLength; ++row){
    for(int column = 0; column < m_columnLength; ++column){
      magnetisation += m_spinMatrix(row, column);
      energy -= m_spinMatrix(row, column)
                * (m_spinMatrix(row, periodic(column-1, m_columnLength))
                    + m_spinMatrix(periodic(row-1, m_rowLength), column)
                  );
    }
  }
  energy *= m_couplingParameter;
  return;
}

void IsingModel::flipSpin(int &row, int &column){
  m_spinMatrix(row, column) *= -1;
}

double IsingModel::calculateDeltaEnergy(int &row, int &column){
  m_deltaEnergy = 2 * m_couplingParameter * m_spinMatrix(row, column)
                * (m_spinMatrix(row, periodic(column-1, m_columnLength))
                    + m_spinMatrix(periodic(row-1, m_rowLength), column)
                    + m_spinMatrix(row, periodic(column+1, m_columnLength))
                    + m_spinMatrix(periodic(row+1, m_rowLength), column)
                  );
  return m_deltaEnergy;
}

double IsingModel::calculateDeltaMagnetisation(int &row, int &column){
  m_deltaMagnetisation = 2.0*m_spinMatrix(row, column);
  return m_deltaMagnetisation;
}
