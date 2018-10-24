#ifndef SYSTEM_HPP
#define	SYSTEM_HPP

#include <vector>
#include <string>
#include <cmath>
#include <fstream>
#include <boost/filesystem.hpp>

#include "planet.hpp"

using namespace std;

// klasse for å løse newtons gravitasjonslov for varierende planetsystemer
class System{
private:
  double m_pi = 3.14159265358979;
  double m_g = 4.0*m_pi*m_pi;
  double m_c = 63239.7263;
  string m_name, m_directory;
  int m_numberofPlanets;
  int m_integrationSteps = 0;
  vector<Planet> m_planets;              // vektor med planetene
  double m_beta = 2;
  bool m_write = false;                  // true hvis man skriver til fil
  vector<ofstream> m_files;              // vektor med filer som åpnes i solver
  bool m_relativistic = false;           // true hvis man regner relativistisk
  bool m_writePerihelion = false;        // sjekker om vi vil gjøre perihelion oppg
  double m_writeParameter = 100;         // skrive ut til fil for hver 100 verdi


  void initPlanets();
  void initFiles();
  void closeFiles();
  Coordinate calculateAcc(int i, int j);

public:
  System()
  : m_name("Sol"), m_planets(vector<Planet>()), m_numberofPlanets(0) {
    // standard bygger
  }
  System(string name, vector<Planet>& planets)
  : m_name(name), m_planets(planets), m_numberofPlanets(planets.size()) {
      // bygger for å initialisere alle kameratparametrene
    }
  // offentlige funksjoner som kan bli brukt utenfor klassen
  string getName(){return m_name;}
  int getNumberofPlanets(){return m_numberofPlanets;}
  void calculateCenterofMass();
  void solveForwardEuler(double endtime, double dt);
  void solveVelocityVerlet(double endtime, double dt);
  void writetoFile(string folder);
  void relativistic();
  void setBeta(double beta){m_beta=beta;}
  void scalePlanetInitVel(double scale, int planet);
  void scalePlanetMass(double scale, int planet);
  double getPlanetInitVel(int planet){return m_planets[planet].m_initVel.norm();}
  double getEnergyChange();
  double getAngularMomentumChange();
  void writePerihelion(string folder);
  void writeForEach(double writeParameter){m_writeParameter = writeParameter;}
};

#endif /* SYSTEM_HPP */
