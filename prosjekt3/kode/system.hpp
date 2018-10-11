#ifndef SYSTEM_HPP
#define	SYSTEM_HPP

#include <string>
#include <cmath>
#include <fstream>

#include "planet.hpp"

using namespace std;

// klasse for å løse newtons gravitasjonslov for varierende planetsystemer
class System{
  // kamerat deklareringer
private:
  double m_pi = 3.14159265358979;
  double m_sunM = 2e+30;
  double m_g = 4.0*m_pi*m_pi/m_sunM;
  string m_name;
  int m_numberofPlanets;
  int m_integrationSteps = 0;
  Planet *m_planets;

  void initPlanets();
  const Planet& copy_planet() const {return *m_planets;}

public:
  System(string name, Planet *planets, int numberofPlanets)
  : m_name(name), m_planets(planets), m_numberofPlanets(numberofPlanets) {
      // bygger for å initialisere alle kameratparametrene
    }
  System(const System& rocks)
  : m_planets(new Planet(rocks.copy_planet())){
    // kopi-bygger for å håndtere flere instanser av klassen
  }
  ~System(){delete[] m_planets;} // ødelegger for å deallokere minne
  // offentlige funksjoner som kan bli brukt utenfor klassen

  string getName(){return m_name;}
  int getNumberofPlanets(){return m_numberofPlanets;}
  void solve(double endtime, double dt, int method, int relativistic);
  Coordinate calculateAcc(int i, int j, int relativistic);
  void velocityVerlet(int i, int j, double dt, int relativistic);
  void forwardEuler(int i, int j, double dt, int relativistic);
  void writetoFile();
};

#endif /* SYSTEM_HPP */
