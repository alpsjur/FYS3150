#include "extractData.hpp"

Planet extract(string filename, int i){
  ifstream inf(filename);
  for (int j=0; j<i; j++){
    string strInput;
    getline(inf, strInput);
  }
  string name, mass, x, y, z, vx, vy, vz;
  inf >> name; inf >> mass;
  inf >> x; inf >> y; inf >> z,
  inf >> vx; inf >> vy; inf >> vz;

  Coordinate initPos(stod(x),stod(y),stod(z));
  Coordinate initVel(stod(vx),stod(vy),stod(vz));
  Planet planet(name, stod(mass)/2e30, initPos, initVel*365.25);

  return planet;
}
