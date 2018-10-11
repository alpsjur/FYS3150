#include "system.hpp"


void System::initPlanets(){
  for(int planet = 0; planet < m_numberofPlanets; ++planet){
    m_planets[planet].initArrays(m_integrationSteps);
    m_planets[planet].initPos();
    m_planets[planet].initVel();
  }
}

//method = 0 gor forward euler, method = 1 gir velocity velocityVerlet
//relativistic = 0 gir klassisk tilfellem relativistic = 1 gir relativistisk tilfelle
void System::solve(double endtime, double dt, int method, int relativistic){
  m_integrationSteps = int(endtime/dt);
  initPlanets();

  for (int i = 0; i < m_integrationSteps-1; ++i){
    //itererer over planetene
    for (int j = 0; j < m_numberofPlanets; ++j){
      if (method == 0){
        forwardEuler(i, j, dt, relativistic);
      }
      if (method == 1){
        velocityVerlet(i, j, dt, relativistic);
      }
    }
  }
}

Coordinate System::calculateAcc(int i, int j,  int relativistic){
  Coordinate force(0,0,0);
  for (int k = 0; k < m_numberofPlanets; ++k){
    if (j != k){
      Coordinate rjk = m_planets[k].m_pos[i] - m_planets[j].m_pos[i];
      if (relativistic == 0){
        force = force + m_g*m_planets[j].m_mass*m_planets[k].m_mass*rjk /
                pow(rjk.norm(),3);
      }
      if (relativistic == 1){
        // kode for relatisvistisk tilfelle
      }
    }
  }
  return force/m_planets[j].m_mass;
}

void System::velocityVerlet(int i, int j, double dt, int relativistic){
  double dt2 = dt/2.0;
  double dtdt2 = dt*dt/2.0;
  Coordinate acc, accNew;

  acc = calculateAcc(i, j, relativistic);
  m_planets[j].m_pos[i+1] = m_planets[j].m_pos[i] + dt*m_planets[j].m_vel[i] + dtdt2*acc;
  accNew = calculateAcc(i+1, j, relativistic);
  m_planets[j].m_vel[i+1] = m_planets[j].m_vel[i] + dt2*(acc + accNew);
}

void System::forwardEuler(int i, int j, double dt, int relativistic){
  Coordinate acc = calculateAcc(i, j, relativistic);
  m_planets[j].m_pos[i+1] = m_planets[j].m_pos[i] + m_planets[j].m_vel[i]*dt;
  m_planets[j].m_vel[i+1] = m_planets[j].m_vel[i] + acc*dt;
}

void System::writetoFile(){
  for (int i = 0; i < m_numberofPlanets; i++){
    string name = m_planets[i].getName();
    ofstream file;
    file.open("../data/" + name + ".dat");
    for (int j = 0; j < m_integrationSteps; j++){
      file << m_planets[i].m_pos[j] << " " << m_planets[i].m_vel[j] << endl;
    }
    file.close();

  }
}
