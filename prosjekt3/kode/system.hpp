#ifndef SYSTEM_HPP
#define	SYSTEM_HPP

#include <string>
#include <cmath>
#include <fstream>
#include <boost/filesystem.hpp>

#include "planet.hpp"

using namespace std;

// klasse for å løse newtons gravitasjonslov for varierende planetsystemer
class System{
  // kamerat deklareringer
private:
  double m_pi = 3.14159265358979;
  double m_g = 4.0*m_pi*m_pi;
  double m_c = 63239.7263;
  string m_name, m_directory;
  int m_numberofPlanets;
  int m_integrationSteps = 0;
  Planet *m_planets;
  double m_beta = 2;
  bool m_write = false;
  ofstream *m_files;
  bool m_relativistic = false;


  void initPlanets();
  const Planet& copyPlanet() const {return *m_planets;}
  void openFiles();
  void closeFiles();

public:
  System(string name, Planet *planets, int numberofPlanets)
  : m_name(name), m_planets(planets), m_numberofPlanets(numberofPlanets) {
      // bygger for å initialisere alle kameratparametrene
    }
  System(const System& rocks)
  : m_planets(new Planet(rocks.copyPlanet())){
    // kopi-bygger for å håndtere flere instanser av klassen
  }
  ~System(){delete[] m_planets, m_files;} // ødelegger for å deallokere minne
  // offentlige funksjoner som kan bli brukt utenfor klassen

  string getName(){return m_name;}
  int getNumberofPlanets(){return m_numberofPlanets;}
  void calculateCenterofMass();
  void solveForwardEuler(double endtime, double dt);
  void solveVelocityVerlet(double endtime, double dt);
  Coordinate calculateAcc(int i, int j);
  void writetoFile(string folder);
  void relativistic();
  void setBeta(double beta){m_beta=beta;}
  void scalePlanetInitVel(double scale, int planet);
  void scalePlanetMass(double scale, int planet);
  double getPlanetInitVel(int planet){return m_planets[planet].m_initVel.norm();}
  double getEnergy();
  double getMomentum();
};

#endif /* SYSTEM_HPP */
