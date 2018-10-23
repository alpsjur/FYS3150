#include "coordinate.hpp"
#include "planet.hpp"

void Planet::initArrays(int integrationSteps){
  m_pos.resize(integrationSteps);
  m_vel.resize(integrationSteps);
}
