#include <fstream>
#include <cstdlib>
#include <iostream>
#include <cmath>

#include "coordinate.hpp"

using namespace std;


class Body {
private:
  string m_name;
  double m_mass;
  Coordinate m_initPos;
  Coordinate m_initVel;
  //Coordinate *m_pos;
  //Coordinate *m_vel;

public:
  Body() : m_name("Navnløs"), m_mass(0), m_initPos(Coordinate()), m_initVel(Coordinate()){
  }
  // A specific constructor
  Body(string name, double mass, Coordinate initPos, Coordinate initVel)
    : m_name(name), m_mass(mass), m_initPos(initPos), m_initVel(initVel){
  }
<<<<<<< HEAD
=======
  Body() : m_name("Navnløs"), m_mass(0), m_initPos(Coordinate()), m_initVel(Coordinate()){
  }
  /*
>>>>>>> af42641f257516f8c4c5635914119453f674d874
  void initArray(int arrayLength){
    m_pos = new Coordinate[arrayLength];
    m_vel = new Coordinate[arrayLength];
  }
  void setPos(int i, Coordinate newPos){
    m_pos[i] = newPos;
  }
  void setVel(int i, Coordinate newVel){
    m_vel[i] = newVel;
  }
  ~Body(){
    delete[] m_pos, m_vel;
  }
  */
};

class System {
private:
  Body m_bodies[];
  int m_numberofBodies;
  int m_arrayLength;
  double m_endtime;
  double m_dt;

public:
<<<<<<< HEAD
  System(Body *bodies, int numberofBodies, double endtime, double dt)
    : m_bodies(bodies), m_numberofBodies(numberofBodies), m_endtime(endtime), m_dt(dt){
    m_arrayLength = (int) m_endtime/m_dt;
    m_bodies = new Body[];
    for (int i = 0; i < m_numberofBodies; ++i){
      //m_bodies[i] = bodies[i];
      m_bodies[i].initArray(m_arrayLength);
    }
  }
  ~System(){
    delete[] m_bodies;
  }
=======
  System(Body bodies[], int numberofBodies, double endtime, double dt)
    :m_numberofBodies(numberofBodies), m_endtime(endtime), m_dt(dt){
    m_arrayLength = (int) m_endtime/m_dt;
    for (int i=0;i<m_numberofBodies;i++){
      m_bodies[i]=bodies[i];
    }
  }
>>>>>>> af42641f257516f8c4c5635914119453f674d874
};

int main(){
  Body jorda("Jorda", 300, Coordinate(1,1,1), Coordinate(2,2,2));
  Body jupiter("Jupiter", 600, Coordinate(2,2,2), Coordinate(1.5,1.5,1.5));
<<<<<<< HEAD

  //Body bodies[] = {jorda, jupiter};
  //System solar(bodies, 2, 10, 0.001);
=======
  Body bodies1[2] = {jorda, jupiter};
  System jordajupiter(bodies1, 2, 10, 0.001);
>>>>>>> af42641f257516f8c4c5635914119453f674d874
  return 0;
}
