#include <string>
#include <cstdlib>

#include "coordinate.hpp"

using namespace std;

class Body{
  // constant/variable declarations
private:
  string m_name;
  const double m_mass;

  Coordinate m_initPos;
  Coordinate m_initVel;

  Coordinate *m_pos;
  Coordinate *m_vel;

public:
  Body(string name, double mass, Coordinate initPos, Coordinate initVel, Coordinate *pos, Coordinate *vel)
    : m_name(name), m_mass(mass), m_initPos(initPos), m_initVel(initVel), m_pos(pos), m_vel(vel){
    // initialising mass, initial position and initial velocity of planet
  }
  // function declarations
  string getName(){return m_name;}
  double getMass(){return m_mass;}
};
