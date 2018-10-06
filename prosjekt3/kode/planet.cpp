
#include "coordinate.hpp"
#include "planet.hpp"

void Planet::initArrays(int integrationSteps){
  m_pos = new Coordinate[integrationSteps];
  m_vel = new Coordinate[integrationSteps];
}
