#include <cstdlib>
#include <string>
#include <iostream>
#include <cmath>

using namespace std;


class Coordinate{
private:
  double m_x;
  double m_y;
  double m_z;
public:
  Coordinate()
  : m_x(0), m_y(0), m_z(0){
    // initialisation if no values are provided
  }
  Coordinate(double x, double y, double z)
  : m_x(x), m_y(y), m_z(z){
    // initialising coordinates if they are provided
  }
  friend bool operator == (const Coordinate &lhs, const Coordinate &rhs);
  friend Coordinate operator + (const Coordinate &lhs, const Coordinate &rhs);
  friend Coordinate operator - (const Coordinate &lhs, const Coordinate &rhs);
  friend Coordinate operator * (double scalar, const Coordinate &rhs);
  friend Coordinate operator * (const Coordinate &lhs, double scalar);
  friend Coordinate operator / (const Coordinate &lhs, double scalar);
  friend ostream &operator << (ostream &out, const Coordinate &rhs);

  double norm(){return sqrt(m_x*m_x + m_y*m_y + m_z*m_z);}
};
