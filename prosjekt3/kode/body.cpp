#include <cstdlib>
#include <iostream>
#include <string>

#include "coordinate.hpp"
#include "body.hpp"

using namespace std;


int main(){
  string nameEarth = "Earth";
  double massEarth = 5.972e24;

  Coordinate initPosEarth(0, 0, 0);
  Coordinate initVelEarth(1, 1, 1);


  Body earth(nameEarth, massEarth, initPosEarth, initVelEarth);
  cout << earth.getMass() << endl;

  string nameJupiter = "Jupiter";
  double massJupiter = 1.898e27;

  Coordinate initPosJupiter(1, 1, 1);
  Coordinate initVelJupiter(2, 2, 2);

  Body jupiter(nameJupiter, massJupiter, initPosJupiter, initVelJupiter);
  cout << jupiter.getMass() << endl;

  return 0;
}
