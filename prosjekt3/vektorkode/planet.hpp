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

  vector<Coordinate> m_pos;
  vector<Coordinate> m_vel;

  // funksjoner som klassen System skal ha tilgang til, men ikke utenfor
  void initArrays(int integrationSteps);
  void initPos(){m_pos[0] = m_initPos;}
  void initVel(){m_vel[0] = m_initVel;}

public:
  Planet()
  : m_name("Planet X"), m_mass(1.0), m_initPos(Coordinate()), m_initVel(Coordinate()) {
    // standard initialisering
  }
  Planet(string name, double mass, Coordinate initPos, Coordinate initVel)
    : m_name(name), m_mass(mass), m_initPos(initPos), m_initVel(initVel){
    // initialiserer alle kameratene og legger vekk minne for pos og vel
  }
  // offentlige funksjoner for bruk utenfor klassen
  Planet& operator = (const Planet& other){
    if(this == &other){return *this;}
    m_name = other.m_name;
    m_mass = other.m_mass;
    m_initPos = other.m_initPos;
    m_initVel = other.m_initVel;
    return *this;
  }
  string getName(){return m_name;}
  double getMass(){return m_mass;}
};

#endif /* PLANET_HPP */
