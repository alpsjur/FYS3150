#ifndef COORDINATE_HPP
#define	COORDINATE_HPP

#include <cstdlib>
#include <string>
#include <iostream>
#include <cmath>

using namespace std;


class Coordinate{
  // klasse for å håndtere vektor aritmetikk. += -= operatorene osv. er ikke implementert enda!
  friend class Planet;
  // kamerat deklareringer
private:
  double m_x;
  double m_y;
  double m_z;
public:
  Coordinate()
  : m_x(0), m_y(0), m_z(0){
    // standard byggeren hvis ingenting blir gitt i initialiseringen. Gir origo
  }
  Coordinate(double x, double y, double z)
  : m_x(x), m_y(y), m_z(z){
    // initialiserer verdiene hvis de blir gitt
  }
  // operator funksjoner for vektor aritmetikk
  friend bool operator == (const Coordinate &lhs, const Coordinate &rhs);
  friend Coordinate operator + (const Coordinate &lhs, const Coordinate &rhs);
  friend Coordinate operator - (const Coordinate &lhs, const Coordinate &rhs);
  friend Coordinate operator * (double scalar, const Coordinate &rhs);
  friend Coordinate operator * (const Coordinate &lhs, double scalar);
  friend double operator * (const Coordinate &lhs, const Coordinate &rhs);
  friend Coordinate operator ^ (const Coordinate &lhs, const Coordinate &rhs);
  friend Coordinate operator / (const Coordinate &lhs, double scalar);
  friend ostream &operator << (ostream &out, const Coordinate &rhs);

  double norm(){return sqrt(m_x*m_x + m_y*m_y + m_z*m_z);}
};

#endif /* COORDINATE_HPP */
