#include <string>

#include "catch.hpp"
#include "coordinate.hpp"
#include "planet.hpp"
#include "system.hpp"
#include "extractData.hpp"


TEST_CASE("ENERGY AND MOMENTUM CONSERVATION"){
  double endtime = 100;
  double dt = 0.001;

  string filename = "../data/body051018.dat";

  string nameSun = "Sun";
  double massSun = 1;
  Coordinate initPosSun(0, 0, 0);
  Coordinate initVelSun(0, 0, 0);
  Planet sun(nameSun, massSun, initPosSun, initVelSun);

  Planet earth = extract(filename, 2);

  vector<Planet> sunEarthList(2);
  sunEarthList[0] = sun;
  sunEarthList[1] = earth;

  System sunEarth("Sun-Earth system", sunEarthList);
  sunEarth.calculateCenterofMass();
  sunEarth.solveVelocityVerlet(endtime, dt);
  double E = sunEarth.getEnergyChange();
  double L = sunEarth.getAngularMomentumChange();

  REQUIRE(E == Approx(0));
  REQUIRE(L == Approx(0));

}
