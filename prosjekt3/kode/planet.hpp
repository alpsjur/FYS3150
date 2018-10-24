#ifndef PLANET_HPP
#define	PLANET_HPP

#include <string>
#include <vector>

#include "coordinate.hpp"

using namespace std;



class Planet{
  friend class System;
  // kamerat deklareringer
private:
  string m_name;
  double m_mass;

  Coordinate m_initPos;
  Coordinate m_initVel;

  Coordinate m_pos;
  Coordinate m_vel;

  // funksjoner som klassen System skal ha tilgang til, men ikke utenfor
  void initPos(){m_pos = m_initPos;}
  void initVel(){m_vel = m_initVel;}

public:
  Planet()
  : m_name("Planet X"), m_mass(1.0), m_initPos(Coordinate()), m_initVel(Coordinate()) {
    // standard initialisering
  }
  Planet(string name, double mass, Coordinate initPos, Coordinate initVel)
    : m_name(name), m_mass(mass), m_initPos(initPos), m_initVel(initVel){
    // initialiserer alle kameratene
  }
  string getName(){return m_name;}
  double getMass(){return m_mass;}
};

#endif /* PLANET_HPP */
