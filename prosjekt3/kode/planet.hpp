#ifndef PLANET_HPP
#define	PLANET_HPP

#include <string>

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

  Coordinate *m_pos;
  Coordinate *m_vel;

  // funksjoner som klassen System skal ha tilgang til, men ikke utenfor
  const Coordinate& copy_coords(Coordinate *coords) const {return *coords;}
  void initArrays(int integrationSteps);
  void initPos(){m_pos[0] = m_initPos;}
  void initVel(){m_vel[0] = m_initVel;}
  void setPos(int i, Coordinate pos){m_pos[i] = pos;}
  void setVel(int i, Coordinate vel){m_vel[i] = vel;}

public:
  Planet()
  : m_name("Planet X"), m_mass(1.0), m_initPos(Coordinate()), m_initVel(Coordinate()){
    // standard initialisering
  }
  Planet(string name, double mass, Coordinate initPos, Coordinate initVel)
    : m_name(name), m_mass(mass), m_initPos(initPos), m_initVel(initVel),
      m_pos(new Coordinate), m_vel(new Coordinate) {
    // initialiserer alle kameratene og legger vekk minne for pos og vel
  }
  Planet(const Planet& other)
  : m_pos(new Coordinate(other.copy_coords(m_pos))),
    m_vel(new Coordinate(other.copy_coords(m_vel))){
    // kopi-bygger for å håndtere flere planeter om gangen
  }
  ~Planet(){delete[] m_pos, m_vel;}  // ødelegger for å deallokere minne
  // offentlige funksjoner for bruk utenfor klassen
  Planet& operator = (const Planet &other){
    if(this == &other){return *this;}
    m_name = other.m_name;
    m_mass = other.m_mass;
    m_initPos = other.m_initPos;
    m_initVel = other.m_initVel;

    return *this;
  }
  string getName(){return m_name;}
  double getMass(){return m_mass;}
  Coordinate getPos(int i){return m_pos[i];}
  Coordinate getVel(int i){return m_vel[i];}
};

#endif /* PLANET_HPP */
