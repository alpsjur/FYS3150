#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include "coordinate.hpp"
#include "planet.hpp"
#include "system.hpp"
#include "extractData.hpp"

using namespace std;

int main(){
  // EKSEMPLER PÅ HVORDAN DET FUNKER
  string filename = "../data/body051018.dat";

  string nameSun = "Sun";
  double massSun = 2e30;

  Coordinate initPosSun(0, 0, 0);
  Coordinate initVelSun(0, 0, 0);

  Planet sun(nameSun, massSun, initPosSun, initVelSun);

  string sysname = "testsystem";
  Planet *planetlist = new Planet[9];
  for (int i = 0; i < 8; i++){
    planetlist[i] = extract(filename, i);
  }
  planetlist[8] = sun;

  System solarsys(sysname, planetlist, 9);
  cout << solarsys.getName() << endl;
  cout << solarsys.getNumberofPlanets()<< endl;
  double endtime = 100;
  double dt = 0.001;
  solarsys.solve(endtime, dt, 1, 0);
  solarsys.writetoFile();
  return 0;
}
