#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include "coordinate.hpp"
#include "planet.hpp"
#include "system.hpp"

using namespace std;

int main(){
  // EKSEMPLER PÅ HVORDAN DET FUNKER
  Coordinate pos(1, 1, 1);
  Coordinate vel(1, 1, 1);
  Coordinate acc(1, 1, 1);
  pos = pos + vel + acc;
  cout << pos << endl;

  string nameEarth = "Earth";
  double massEarth = 5.972e24;

  Coordinate initPosEarth(0, 0, 0);
  Coordinate initVelEarth(1, 1, 1);


  Planet earth(nameEarth, massEarth, initPosEarth, initVelEarth);
  cout << earth.getMass() << endl;

  string nameJupiter = "Jupiter";
  double massJupiter = 1.898e27;

  Coordinate initPosJupiter(1, 1, 1);
  Coordinate initVelJupiter(2, 2, 2);

  Planet jupiter(nameJupiter, massJupiter, initPosJupiter, initVelJupiter);
  cout << jupiter.getMass() << endl;


  string sysname = "testsystem";
  Planet *planetlist = new Planet[2];
  planetlist[0] = earth;
  planetlist[1] = jupiter;

  System earthjupiter(sysname, planetlist, 2);
  cout << earthjupiter.getName() << endl;
  cout << earthjupiter.getNumberofPlanets()<< endl;
  double endtime = 10;
  double dt = 0.001;
  earthjupiter.solve(endtime, dt);

  return 0;
}
