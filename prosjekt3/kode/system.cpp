#include "system.hpp"


void System::initPlanets(){
  for(int planet = 0; planet < m_numberofPlanets; ++planet){
    m_planets[planet].initArrays(m_integrationSteps);
    m_planets[planet].initPos();
    m_planets[planet].initVel();
  }
}

void System::calculateCenterofMass(){
  double totMass = 0;
  Coordinate sum(0,0,0);
  Coordinate planetMomentum(0,0,0);
  for (int i = 0; i < m_numberofPlanets ; i++){
    totMass = totMass + m_planets[i].m_mass;
    sum = sum + m_planets[i].m_mass * m_planets[i].m_initPos;
    planetMomentum = planetMomentum + m_planets[i].m_mass * m_planets[i].m_initVel;
  }
  Coordinate centrum = sum/totMass;
  m_planets[0].m_initVel = m_planets[0].m_initVel-planetMomentum;
  for (int i = 0; i < m_numberofPlanets; i++){
    m_planets[i].m_initPos = m_planets[i].m_initPos - centrum;
  }

}


void System::solveEuler(double endtime, double dt){
  m_integrationSteps = int(endtime/dt);
  initPlanets();
  for (int i = 0; i < m_integrationSteps-1; ++i){
    //itererer over planetene
    for (int j = 0; j < m_numberofPlanets; ++j){
        forwardEuler(i, j, dt);
    }
  }
}

void System::solveVerlet(double endtime, double dt){
  m_integrationSteps = int(endtime/dt);
  initPlanets();

  for (int i = 0; i < m_integrationSteps-1; ++i){
    //itererer over planetene
    for (int j = 0; j < m_numberofPlanets; ++j){
        velocityVerlet(i, j, dt);
    }
  }
}

Coordinate System::calculateAcc(int i, int j){
  Coordinate force(0,0,0);
  Coordinate forcejk(0,0,0);
  double correction = 1;
  for (int k = 0; k < m_numberofPlanets; ++k){
    if (j != k){
      Coordinate rjk = m_planets[k].m_pos[i] - m_planets[j].m_pos[i];
      forcejk=  m_g*m_planets[j].m_mass*m_planets[k].m_mass*rjk /
                pow(rjk.norm(),m_beta+1);
      if (m_relativistic == 1){
        double l = (rjk^m_planets[j].m_vel[i]).norm();
        double r = rjk.norm();
        correction = 1 + 3*pow(l,2)/pow(r*m_c,2);
      }
      force = force + forcejk*correction;
    }
  }
  return force/m_planets[j].m_mass;
}

void System::velocityVerlet(int i, int j, double dt){
  double dt2 = dt/2.0;
  double dtdt2 = dt*dt/2.0;
  Coordinate acc, accNew;

  acc = calculateAcc(i, j);
  m_planets[j].m_pos[i+1] = m_planets[j].m_pos[i] + dt*m_planets[j].m_vel[i] + dtdt2*acc;
  accNew = calculateAcc(i+1, j);
  m_planets[j].m_vel[i+1] = m_planets[j].m_vel[i] + dt2*(acc + accNew);
}

void System::forwardEuler(int i, int j, double dt){
  Coordinate acc = calculateAcc(i, j);
  m_planets[j].m_pos[i+1] = m_planets[j].m_pos[i] + m_planets[j].m_vel[i]*dt;
  m_planets[j].m_vel[i+1] = m_planets[j].m_vel[i] + acc*dt;
}

void System::writetoFile(string folder){
  for (int i = 0; i < m_numberofPlanets; i++){
    string name = m_planets[i].getName();
    ofstream file;
    file.open(folder+ "/" + name + ".dat");
    for (int j = 0; j < m_integrationSteps; j++){
      file << m_planets[i].m_pos[j] << " " << m_planets[i].m_vel[j] << endl;
    }
    file.close();

  }
}

void System::relativistic(string arg){
  if(arg == "on"){
    m_relativistic=1;
  }
  if(arg=="off"){
    m_relativistic=0;
  }
}

void System::scalePlanetInitVel(double scale, int planet){
  m_planets[planet].m_initVel=m_planets[planet].m_initVel*scale;
}

void System::scalePlanetMass(double scale, int planet){
  m_planets[planet].m_mass = m_planets[planet].m_mass*scale;
}
