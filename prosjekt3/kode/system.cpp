
#include "system.hpp"

void System::initPlanets(int integrationSteps){
  for(int planet = 0; planet < m_numberofPlanets; ++planet){
    m_planets[planet].initArrays(integrationSteps);
    m_planets[planet].initPos();
    m_planets[planet].initVel();
  }
}

void System::solve(double endtime, double dt){
  int integrationSteps = int(endtime/dt);
  initPlanets(integrationSteps);
  // KODE FOR Å LØSE NEWTONS LOVER
}
