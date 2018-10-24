#include "system.hpp"


void System::initPlanets(){
  for(int planet = 0; planet < m_numberofPlanets; ++planet){
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
  double local_writeParameter;            // denne sparer noen flops

  Coordinate acc;

  initPlanets();
  if(m_write){
    initFiles();
    // itererer over tidssteg
    for (int i = 0; i < m_integrationSteps-1; ++i){
      //itererer over planetene
      local_writeParameter = (i+1)/m_writeParameter;
      for (int j = 0; j < m_numberofPlanets; ++j){
        acc = calculateAcc(i, j);
        m_planets[j].m_pos = m_planets[j].m_pos + m_planets[j].m_vel*dt;
        m_planets[j].m_vel = m_planets[j].m_vel + acc*dt;
        if(local_writeParameter == int(local_writeParameter)){
          m_files[j] << m_planets[j].m_pos << " " << m_planets[j].m_vel << endl;
        }
      }
    }
    closeFiles();
  }
  else{
    for (int i = 0; i < m_integrationSteps-1; ++i){
      //itererer over planetene
      for (int j = 0; j < m_numberofPlanets; ++j){
        acc = calculateAcc(i, j);
        m_planets[j].m_pos = m_planets[j].m_pos + m_planets[j].m_vel*dt;
        m_planets[j].m_vel = m_planets[j].m_vel + acc*dt;
      }
    }
  }
}

void System::solveVelocityVerlet(double endtime, double dt){
  m_integrationSteps = int(endtime/dt);
  double dt2 = dt/2.0;
  double dtdt2 = dt*dt/2.0;
  double local_writeParameter;       // denne blir brukt til å spare flops

  vector<Coordinate> acc(m_numberofPlanets);
  Coordinate accNew;
  initPlanets();
  if(m_write){
    initFiles();
    // itererer over tid
    for (int i = 0; i < m_integrationSteps-1; ++i){
      //itererer over planetene
      local_writeParameter = (i+1)/m_writeParameter;
      for (int j = 0; j < m_numberofPlanets; ++j){
        acc[j] = calculateAcc(i, j);
        m_planets[j].m_pos = m_planets[j].m_pos + dt*m_planets[j].m_vel + dtdt2*acc[j];
      }
      for (int j = 0; j < m_numberofPlanets; ++j){
        accNew = calculateAcc(i+1, j);
        m_planets[j].m_vel = m_planets[j].m_vel + dt2*(acc[j] + accNew);
        if(local_writeParameter == int(local_writeParameter)){
          m_files[j] << m_planets[j].m_pos << " " << m_planets[j].m_vel << endl;
        }
      }
    }
    closeFiles();
  }
  // her var vi litt desp, burde ordne dette i python istedet
  else if(m_writePerihelion){
    double r;
    vector<double> rp(m_numberofPlanets-1);
    vector<double> rpp(m_numberofPlanets-1);
    initFiles();
    for (int i = 0; i < m_integrationSteps-1; ++i){
      //itererer over planetene
      for (int j = 0; j < m_numberofPlanets; ++j){
        if(j!=0){
          Coordinate r0j = m_planets[j].m_pos - m_planets[0].m_pos;
          r = r0j.norm();
          if((r>rp[j-1]) && (rpp[j-1]>rp[j-1])){
            m_files[j-1] << m_planets[j].m_pos << " " << i << endl;
          }
          rpp[j-1] = rp[j-1]; rp[j-1] = r;
        }
        acc[j] = calculateAcc(i, j);
        m_planets[j].m_pos = m_planets[j].m_pos + dt*m_planets[j].m_vel + dtdt2*acc[j];
      }
      for (int j = 0; j < m_numberofPlanets; ++j){
        accNew = calculateAcc(i+1, j);
        m_planets[j].m_vel = m_planets[j].m_vel + dt2*(acc[j] + accNew);
      }
    }
    closeFiles();
  }

  else{
    for (int i = 0; i < m_integrationSteps-1; ++i){
      //itererer over planetene
      for (int j = 0; j < m_numberofPlanets; ++j){
        acc[j] = calculateAcc(i, j);
        m_planets[j].m_pos = m_planets[j].m_pos + dt*m_planets[j].m_vel + dtdt2*acc[j];
      }
      for (int j = 0; j < m_numberofPlanets; ++j){
        accNew = calculateAcc(i+1, j);
        m_planets[j].m_vel = m_planets[j].m_vel + dt2*(acc[j] + accNew);
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
      Coordinate rjk = m_planets[k].m_pos - m_planets[j].m_pos;
      forcejk =  m_g*m_planets[j].m_mass*m_planets[k].m_mass*rjk/pow(rjk.norm(),m_beta+1);
      if (m_relativistic){
        double l = (rjk^m_planets[j].m_vel).norm();
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
  if(m_writePerihelion){
    m_files.resize(m_numberofPlanets-1);
    for (int j = 0; j < m_numberofPlanets-1; ++j){
      string filename;
      filename = m_directory + "/" + m_planets[j+1].getName() + "Perihelion.dat";
      m_files[j].open(filename);
    }
  }
  else{
    m_files.resize(m_numberofPlanets);
    for (int j = 0; j < m_numberofPlanets; ++j){
      string filename;
      filename = m_directory + "/" + m_planets[j].getName() + ".dat";
      m_files[j].open(filename);
      m_files[j] << m_planets[j].m_initPos << " " << m_planets[j].m_initVel << endl;
    }
  }

}

void System::closeFiles(){
  if(m_writePerihelion){
    for(int j = 0; j < m_numberofPlanets-1; ++j){
      m_files[j].close();
    }
  }
  else{
    for(int j = 0; j < m_numberofPlanets; ++j){
      m_files[j].close();
    }
  }
}

void System::writetoFile(string folder){
  m_write = true;
  m_directory = "../data/" + folder;
}

double System::getEnergyChange(){
  // denne funksjonen er bare korrekt for to legemer akkurat nå
  double dE = 0;
  double E0, E1;
  double U0, U1;
  double T0, T1;
  // vi anser ikke den totale energien til sola
  for(int j = 1; j < m_numberofPlanets; ++j){
    U0 = -m_g*m_planets[j].getMass()/m_planets[j].m_initPos.norm();
    T0 = 0.5*m_planets[j].getMass()*m_planets[j].m_initVel*m_planets[j].m_initVel;
    E0 = U0 + T0;
    U1 = -m_g*m_planets[j].getMass()/m_planets[j].m_pos.norm();
    T1 = 0.5*m_planets[j].getMass()*m_planets[j].m_vel*m_planets[j].m_vel;
    E1 = U1 + T1;
    dE += E1 - E0;
  }
  return dE;
}

double System::getAngularMomentumChange(){
  // denne funksjonen funker også bare for to legemer
  Coordinate dL, L0, L1;
  for(int j = 1; j < m_numberofPlanets; ++j){
    L0 = m_planets[j].m_initPos^m_planets[j].getMass()*m_planets[j].m_initVel;
    L1 = m_planets[j].m_pos^m_planets[j].getMass()*m_planets[j].m_vel;
    dL = dL + (L1 - L0);
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

void System::writePerihelion(string folder){
  m_writePerihelion = true;
  m_write = false;
  m_directory = "../data/" + folder;
}
