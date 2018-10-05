#include <string>
#include <cmath>

using namespace std;

class System{
  // member declarations
private:
  double m_pi = 3.14159265358979;
  double m_gm = 4.0*m_pi*m_pi;
  string m_name;
  int m_numberofBodies;
  Body *m_bodies;

  void initArrays();

public:
  System(string name, Body bodies[])
  : m_name(name), m_numberofBodies(sizeof(bodies)/sizeof(*bodies)),
    m_bodies(bodies) {
      // constructor to initialise variables and constants
    }
  //System(const System& rocks)
  //: m_bodies(new Body(rocks.copy_body())) {
    // copy constructor to handle multiple system instances
  //}
  //~System(){delete[] m_bodies;}  // destructor to deallocate m_bodies
  // public functions
  //const Body& copy_body() const {return *m_bodies;}
  void solve(double endtime, double dt);

};
