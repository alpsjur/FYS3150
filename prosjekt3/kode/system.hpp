#ifndef SYSTEM_HPP
#define	SYSTEM_HPP

#include <string>
#include <cmath>

#include "planet.hpp"

using namespace std;

// klasse for å løse newtons gravitasjonslov for varierende planetsystemer
class System{
  // kamerat deklareringer
private:
  double m_pi = 3.14159265358979;
  double m_gm = 4.0*m_pi*m_pi;
  string m_name;
  int m_numberofPlanets;
  Planet *m_planets;

  void initPlanets(int integrationSteps);
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
  void solve(double endtime, double dt);
};

#endif /* SYSTEM_HPP */
