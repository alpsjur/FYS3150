#include <string>
#include <cstdlib>

using namespace std;

class Body{
  // constant/variable declarations
private:
  string m_name;
  double m_mass;

  Coordinate m_initPos;
  Coordinate m_initVel;

  Coordinate *m_pos;
  Coordinate *m_vel;

public:
  Body()
  : m_name("UFO"), m_mass(1.0), m_initPos(Coordinate()), m_initVel(Coordinate()),
  m_pos(new Coordinate), m_vel(new Coordinate) {
    // default initialisation
  }
  Body(string name, double mass, Coordinate initPos, Coordinate initVel)
    : m_name(name), m_mass(mass), m_initPos(initPos), m_initVel(initVel),
      m_pos(new Coordinate), m_vel(new Coordinate) {
    // initialising mass, initial position and initial velocity of planet
  }
  Body(const Body& rock)
  : m_pos(new Coordinate(rock.copy_coord(m_pos))),
    m_vel(new Coordinate(rock.copy_coord(m_vel))) {
    // copy constructor for handling multiple instances of class
  }
  ~Body(){delete[] m_pos, m_vel;}  // destructor freeing memory for pointers
  // function declarations
  const Coordinate& copy_coord(Coordinate *coord) const {return *coord;}
  string getName(){return m_name;}
  double getMass(){return m_mass;}
  Coordinate getInitPos(){return m_initPos;}
  Coordinate getInitVel(){return m_initVel;}
};
