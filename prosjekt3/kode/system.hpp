#ifndef SYSTEM_HPP
#define	SYSTEM_HPP

#include <string>
#include <cmath>

#include "planet.hpp"

using namespace std;

class System{
  // member declarations
private:
  double m_pi = 3.14159265358979;
  double m_gm = 4.0*m_pi*m_pi;
  string m_name;
  int m_numberofPlanets;
  Planet *m_planets;

public:
  System(string name, Planet *planets, int numberofPlanets)
  : m_name(name), m_planets(planets), m_numberofPlanets(numberofPlanets) {
      // constructor to initialise variables and constants
    }
  System(const System& rocks)
  : m_planets(new Planet(rocks.copy_planet())){
    // copy constructor
  }
  ~System(){delete[] m_planets;}
  // public functions
  const Planet& copy_planet() const {return *m_planets;}

  string getName(){return m_name;}
  int getNumberofPlanets(){return m_numberofPlanets;}
  void initPlanets(int integrationSteps);
  void solve(double endtime, double dt);
};

#endif /* SYSTEM_HPP */
