#include <iostream>
#include <fstream>
#include <string>
#include "coordinate.hpp"
#include "body.hpp"

using namespace std;

Body extract(string filename, int i){
  ifstream inf(filename);
  while(j=0;j<i;j++){
    string strInput;
    getline(inf, strInput);
  }
  string name, mass, x, y, z, vx, vy, vz;
  inf >> name; inf >> mass;
  inf >> x; inf >> y; inf >> z,
  inf >> vx; inf >> vy; inf >> vz;

  double mass = atof(mass);
  Coordinate initPos(atof(x),atof(y),atof(z));
  Coordinate initVel(atof(x),atof(vy),atof(vz));
  Body planet(name, mass, initPos, initVel);

  return planet;
}

int main() {
  /* code */
  return 0;
}
