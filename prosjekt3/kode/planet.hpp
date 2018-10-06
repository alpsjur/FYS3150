#ifndef BODY_HPP
#define	BODY_HPP

#include <string>

#include "coordinate.hpp"

using namespace std;



class Planet{
  // member declarations
private:
  string m_name;
  double m_mass;

  Coordinate m_initPos;
  Coordinate m_initVel;

  Coordinate *m_pos;
  Coordinate *m_vel;


public:
  Planet()
  : m_name("Planet X"), m_mass(1.0), m_initPos(Coordinate()), m_initVel(Coordinate()),
  m_pos(new Coordinate), m_vel(new Coordinate) {
    // default initialisation
  }
  Planet(string name, double mass, Coordinate initPos, Coordinate initVel)
    : m_name(name), m_mass(mass), m_initPos(initPos), m_initVel(initVel),
      m_pos(new Coordinate), m_vel(new Coordinate){
    // initialising mass, initial position and initial velocity of planet
  }
  Planet(const Planet& rock)
  : m_pos(new Coordinate(rock.copy_coord(m_pos))),
    m_vel(new Coordinate(rock.copy_coord(m_vel))){
    // copy constructor for handling multiple instances of class
  }
  ~Planet(){delete[] m_pos, m_vel;}  // destructor freeing memory for pointers
  // function declarations
  const Coordinate& copy_coord(Coordinate *coord) const {return *coord;}
  string getName(){return m_name;}
  double getMass(){return m_mass;}
  //Coordinate getPos(int i){return m_pos[i];}
  //Coordinate getVel(int i){return m_vel[i];}
  //Coordinate setPos(int i, Coordinate pos){m_pos[i] = pos;}
  //Coordinate setVel(int i, Coordinate vel){m_vel[i] = vel;}
  void initArrays(int integrationSteps);
  void initPos(){m_pos[0] = m_initPos;}
  void initVel(){m_vel[0] = m_initVel;}
};

#endif /* BODY_HPP */
