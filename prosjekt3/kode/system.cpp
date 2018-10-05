
#include "coordinate.hpp"
#include "body.hpp"
#include "system.hpp"

void System::initArrays(int integrationSteps){

  for(int body 0; body < m_numberofBodies; ++body){
    m_bodies[body].m_pos = new Coordinate[integrationSteps];
    m_bodies[body].m_vel = new Coordinate[integrationSteps];
    m_bodies[body].m_pos[0] = m_bodies[body].getInitPos();
    m_bodies[body].m_vel[0] = m_bodies[body].getInitVel();
  }
}

void System::solve(double endtime, double dt){
  double integrationSteps = int(endtime/dt);

  return;
}


int main(){

  return 0;
}
