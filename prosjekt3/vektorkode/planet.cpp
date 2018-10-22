#include "coordinate.hpp"
#include "planet.hpp"


void Planet::initArrays(int integrationSteps){
  m_pos.reserve(integrationSteps);
  m_vel.reserve(integrationSteps);
}
