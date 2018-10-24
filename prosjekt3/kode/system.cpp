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
  for (int i = 0; i < m_numberofPlanets ; ++i){
    totMass = totMass + m_planets[i].m_mass;
    sum = sum + m_planets[i].m_mass * m_planets[i].m_initPos;
    planetMomentum = planetMomentum + m_planets[i].m_mass * m_planets[i].m_initVel;
  }
  Coordinate centrum = sum/totMass;
  m_planets[0].m_initVel = m_planets[0].m_initVel-planetMomentum;
  for (int i = 0; i < m_numberofPlanets; ++i){
    m_planets[i].m_initPos = m_planets[i].m_initPos - centrum;
  }
}


void System::solveForwardEuler(double endtime, double dt){
  m_integrationSteps = int(endtime/dt);
  Coordinate acc;

  initPlanets();
  if(m_write){
    initFiles();
    for (int i = 0; i < m_integrationSteps-1; ++i){
      //itererer over planetene
      for (int j = 0; j < m_numberofPlanets; ++j){
        acc = calculateAcc(i, j);
        m_planets[j].m_pos[i+1] = m_planets[j].m_pos[i] + m_planets[j].m_vel[i]*dt;
        m_planets[j].m_vel[i+1] = m_planets[j].m_vel[i] + acc*dt;
        m_files[j] << m_planets[j].m_pos[i+1] << " " << m_planets[j].m_vel[i+1] << endl;
      }
    }
    closeFiles();
  }
  else{
    for (int i = 0; i < m_integrationSteps-1; ++i){
      //itererer over planetene
      for (int j = 0; j < m_numberofPlanets; ++j){
        acc = calculateAcc(i, j);
        m_planets[j].m_pos[i+1] = m_planets[j].m_pos[i] + m_planets[j].m_vel[i]*dt;
        m_planets[j].m_vel[i+1] = m_planets[j].m_vel[i] + acc*dt;
      }
    }
  }
}

void System::solveVelocityVerlet(double endtime, double dt){
  m_integrationSteps = int(endtime/dt);
  double dt2 = dt/2.0;
  double dtdt2 = dt*dt/2.0;
  vector<Coordinate> acc;
  acc.reserve(m_numberofPlanets);
  Coordinate accNew;

  initPlanets();
  if(m_write){
    initFiles();
    for (int i = 0; i < m_integrationSteps-1; ++i){
      //itererer over planetene
      for (int j = 0; j < m_numberofPlanets; ++j){
        acc[j] = calculateAcc(i, j);
        m_planets[j].m_pos[i+1] = m_planets[j].m_pos[i] + dt*m_planets[j].m_vel[i] + dtdt2*acc[j];
      }
      for (int j = 0; j < m_numberofPlanets; ++j){
        accNew = calculateAcc(i+1, j);
        m_planets[j].m_vel[i+1] = m_planets[j].m_vel[i] + dt2*(acc[j] + accNew);
        m_files[j] << m_planets[j].m_pos[i+1] << " " << m_planets[j].m_vel[i+1] << endl;
      }
    }
    closeFiles();
  }
  else{
    for (int i = 0; i < m_integrationSteps-1; ++i){
      //itererer over planetene
      for (int j = 0; j < m_numberofPlanets; ++j){
        acc[j] = calculateAcc(i, j);
        m_planets[j].m_pos[i+1] = m_planets[j].m_pos[i] + dt*m_planets[j].m_vel[i] + dtdt2*acc[j];
      }
      for (int j = 0; j < m_numberofPlanets; ++j){
        accNew = calculateAcc(i+1, j);
        m_planets[j].m_vel[i+1] = m_planets[j].m_vel[i] + dt2*(acc[j] + accNew);
      }
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
      forcejk =  m_g*m_planets[j].m_mass*m_planets[k].m_mass*rjk/pow(rjk.norm(),m_beta+1);
      if (m_relativistic){
        double l = (rjk^m_planets[j].m_vel[i]).norm();
        double r = rjk.norm();
        correction = 1 + 3*l*l/pow(r*m_c,2);
      }
      force = force + forcejk*correction;
    }
  }
  return force/m_planets[j].m_mass;
}

void System::initFiles(){
  if(boost::filesystem::exists(m_directory) == false){
    boost::filesystem::create_directories(m_directory);
  }
  m_files.resize(m_numberofPlanets);
  for (int j = 0; j < m_numberofPlanets; ++j){
    string filename;
    filename = m_directory + "/" + m_planets[j].getName() + ".dat";
    m_files[j].open(filename);
    m_files[j] << m_planets[j].m_pos[0] << " " << m_planets[j].m_vel[0] << endl;
  }
}

void System::closeFiles(){
  for(int j = 0; j < m_numberofPlanets; ++j){
    m_files[j].close();
  }
}

void System::writetoFile(string folder){
  m_write = true;
  m_directory = "../data/" + folder;
}

double System::getEnergyTotal(){
  // denne funksjonen er bare korrekt for to legemer akkurat nå
  double dE = 0;
  double E0, E1;
  double U0, U1;
  double T0, T1;
  // vi anser ikke den totale energien til sola
  for(int j = 1; j < m_numberofPlanets; ++j){
    U0 = -m_g*m_planets[j].getMass()/m_planets[j].m_pos[0].norm();
    T0 = 0.5*m_planets[j].getMass()*m_planets[j].m_vel[0]*m_planets[j].m_vel[0];
    E0 = U0 + T0;
    for(int i = 0; i < m_integrationSteps-1; ++i){
      U1 = -m_g*m_planets[j].getMass()/m_planets[j].m_pos[i+1].norm();
      T1 = 0.5*m_planets[j].getMass()*m_planets[j].m_vel[i+1]*m_planets[j].m_vel[i+1];
      E1 = U1 + T1;
      dE += E1 - E0;
    }
  }
  return dE;
}

double System::getAngularMomentumTotal(){
  // denne funksjonen funker også bare for to legemer
  Coordinate dL, L0, L1;

  for(int j = 1; j < m_numberofPlanets; ++j){
    L0 = m_planets[j].m_pos[0]^m_planets[j].getMass()*m_planets[j].m_vel[0];
    for(int i = 0; i < m_integrationSteps-1; ++i){
      L1 = m_planets[j].m_pos[i+1]^m_planets[j].getMass()*m_planets[j].m_vel[i+1];
      dL = dL + (L1 - L0);
    }
  }
  return dL.norm();
}


void System::relativistic(){
  m_relativistic = true;
}

void System::scalePlanetInitVel(double scale, int planet){
  m_planets[planet].m_initVel = m_planets[planet].m_initVel*scale;
}

void System::scalePlanetMass(double scale, int planet){
  m_planets[planet].m_mass = m_planets[planet].m_mass*scale;
}

void System::calculatePelihelion(){
  double r;
  double rp = 1e10;
  double rpp = 1e10;
  if(m_writePelihelion){
    if(boost::filesystem::exists(m_directory) == false){
      boost::filesystem::create_directories(m_directory);
    }
  }
  ofstream file;
  for (int j=1; j < m_numberofPlanets; j++){
    if(m_writePelihelion){
      string filename;
      filename = m_directory + "/" + m_planets[j].getName() + "Pelihelon.dat";
      file.open(filename);
    }
    for (int i=0; i < m_integrationSteps; i++){
      Coordinate r0j = m_planets[j].m_pos[i] - m_planets[0].m_pos[i];
      r = r0j.norm();
      if((r>rp)&&(rpp>rp)){
        if(m_writePelihelion){
          file <<m_planets[j].m_pos[i-1] << " " << i-1 << endl;
        }
      }
      rpp = rp; rp = r;
    }
    if(m_writePelihelion){
      file.close();
    }
  }
}

void System::writePelihelion(string folder){
  m_writePelihelion = true;
  m_directory = "../data/" + folder;
}
