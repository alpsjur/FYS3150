
#include "coordinate.hpp"
#include "planet.hpp"

void Planet::initArrays(int integrationSteps){
  delete[] m_pos, m_vel;
  m_pos = new Coordinate[integrationSteps];
  m_vel = new Coordinate[integrationSteps];
}
