#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include <ctime>

#include "coordinate.hpp"
#include "planet.hpp"
#include "system.hpp"
#include "extractData.hpp"


using namespace std;

void timeAlgorithms(vector<Planet>&, double, double);
void compareEulerVerlet(vector<Planet>&, double, double);
void varyVelocity(vector<Planet>&, double, double, double);
void varyBeta(vector<Planet>&, double, double, double);
void solveEarthJupiter(vector<Planet>&, double, double);
void solveCenterofMass(vector<Planet>&, double, double);
void solveSolarSystem(vector<Planet>&, double, double);
void solveMercuryPrecession(vector<Planet>&, double, double);

int main(int argc, char* argv[]){
  double scenario = atof(argv[1]);
  double endtime = atof(argv[2]);
  double dt = atof(argv[3]);
  //data hentet fra NASA, posisjon og hastighet 05.10.18
  string filename = "../data/body051018.dat";

  //intisialiserer Sola og planetene
  string nameSun = "Sun";
  double massSun = 1;
  Coordinate initPosSun(0, 0, 0);
  Coordinate initVelSun(0, 0, 0);
  Planet sun(nameSun, massSun, initPosSun, initVelSun);

  string nameMercury = "Mercury";
  double massMercury = 3.3E23/2E30;
  Coordinate initPosMercury(0.3075, 0, 0);
  Coordinate initVelMercury(0, 12.44, 0);
  Planet mercury(nameMercury, massMercury, initPosMercury, initVelMercury);

  Planet earth = extract(filename, 2);
  Planet jupiter = extract(filename, 4);

  //initsialiserer planetlistene
  vector<Planet> allplanets(9);
  allplanets[0] = sun;
  for(int i = 0; i<8; ++i){
    allplanets[i+1] = extract(filename, i);
  }

  vector<Planet> sunEarthList(2);
  sunEarthList[0] = sun;
  sunEarthList[1] = earth;


  vector<Planet> sunEarthJupiterList(3);
  sunEarthJupiterList[0] = sun;
  sunEarthJupiterList[1] = earth;
  sunEarthJupiterList[2] = jupiter;

  vector<Planet> sunMercuryList(2);
  sunMercuryList[0] = sun;
  sunMercuryList[1] = mercury;

  if(scenario == 0){timeAlgorithms(sunEarthList, endtime, dt);}
  if(scenario == 1){compareEulerVerlet(sunEarthList, endtime, dt);}
  if(scenario == 2){
    double velocityscale = atof(argv[4]);
    varyVelocity(sunEarthList, endtime, dt, velocityscale);
  }
  if(scenario == 3){
    double beta = atof(argv[4]);
    varyBeta(sunEarthList, endtime, dt, beta);
  }
  if(scenario == 4){solveEarthJupiter(sunEarthJupiterList, endtime, dt);}
  if(scenario == 5){solveCenterofMass(sunEarthJupiterList, endtime, dt);}
  if(scenario == 6){solveSolarSystem(allplanets, endtime, dt);}
  if(scenario == 7){solveMercuryPrecession(sunMercuryList, endtime, dt);}

  return 0;
}

void timeAlgorithms(vector<Planet>& sunEarthList, double endtime, double dt){

  System sunEarth("Sun-Earth system", sunEarthList);

  clock_t c_start = clock();
  sunEarth.solveForwardEuler(endtime, dt);
  clock_t c_end = clock();

  // Beregner CPU-tid i milisekunder
  double eulerTime = 1000.0 * (c_end - c_start) / CLOCKS_PER_SEC;

  c_start = clock();
  sunEarth.solveVelocityVerlet(endtime, dt);
  c_end = clock();

  // Beregner CPU-tid i milisekunder
  double verletTime = 1000.0 * (c_end - c_start) / CLOCKS_PER_SEC;

  cout << eulerTime << " " << verletTime << endl;
}

void compareEulerVerlet(vector<Planet>& sunEarthList, double endtime, double dt){
  System sunEarth("Sun-Earth system", sunEarthList);

  sunEarth.writetoFile("euler_vs_verlet/euler");
  sunEarth.writeForEach(10);
  sunEarth.solveForwardEuler(endtime, dt);

  sunEarth.writetoFile("euler_vs_verlet/verlet");
  sunEarth.writeForEach(10);
  sunEarth.solveVelocityVerlet(endtime, dt);
}

void varyBeta(vector<Planet>& sunEarthList, double endtime, double dt, double beta){

  System sunEarthScale("Sun-Earth system", sunEarthList);

  sunEarthScale.setBeta(beta);                               //endrer beta i kraftfunksjonen
  sunEarthScale.writetoFile("change_beta");
  sunEarthScale.solveVelocityVerlet(endtime, dt);

}


void varyVelocity(vector<Planet>& sunEarthList, double endtime, double dt, double velocityScale){

  System sunEarthScale("Sun-Earth system", sunEarthList);
  sunEarthScale.writeForEach(10);
  sunEarthScale.scalePlanetInitVel(velocityScale, 1);        //skalerer hastigheten til Jorda
  sunEarthScale.writetoFile("escape_velocity");
  sunEarthScale.solveVelocityVerlet(endtime, dt);
  double scaledVelocity = sunEarthScale.getPlanetInitVel(1); //henter skalert initsialhastighet

  cout << scaledVelocity << endl;
}

void solveEarthJupiter(vector<Planet> &sunEarthJupiterList, double endtime, double dt){

  System sunEarthJupiter("Sun-Earth_jupiter system", sunEarthJupiterList);


  //løser for Jupiters masse skalert med 1
  sunEarthJupiter.writetoFile("sun_earth_jupiter/jupiter_mass_1");
  sunEarthJupiter.solveVelocityVerlet(endtime, dt);

  //løser for jupiters masse skalert med 10
  sunEarthJupiter.scalePlanetMass(10,2);
  sunEarthJupiter.writetoFile("sun_earth_jupiter/jupiter_mass_10");
  sunEarthJupiter.solveVelocityVerlet(endtime, dt);

  //løser for jupiters masse skalert med 1000 = 10*100
  sunEarthJupiter.scalePlanetMass(100,2);
  sunEarthJupiter.writetoFile("sun_earth_jupiter/jupiter_mass_1000");
  sunEarthJupiter.solveVelocityVerlet(endtime, dt);
}

void solveCenterofMass(vector<Planet> &sunEarthJupiterList, double endtime, double dt){
  System sunEarthJupiter("Sun-Earth-jupiter system", sunEarthJupiterList);

  //løser uten å sette massesenteret som origo
  sunEarthJupiter.writetoFile("sun_earth_jupiter/sun_origo");
  sunEarthJupiter.solveVelocityVerlet(endtime, dt);

  //løser etter å ha satt massesenteret som origo og gitt sola initsialhastighet
  sunEarthJupiter.calculateCenterofMass();
  sunEarthJupiter.writetoFile("sun_earth_jupiter/mass_origo");
  sunEarthJupiter.solveVelocityVerlet(endtime, dt);
}

void solveSolarSystem(vector<Planet> &allplanets, double endtime, double dt){

  System solarsystem("Solar system", allplanets);
  solarsystem.calculateCenterofMass();
  solarsystem.writetoFile("solarsystem");
  solarsystem.solveVelocityVerlet(endtime, dt);
}

void solveMercuryPrecession(vector<Planet> &sunMercuryList, double endtime, double dt){
  System sunMercuryClassical("Sun-Mercury classical system", sunMercuryList);
  System sunMercuryRelativistic("Sun-Mercury relativistic system", sunMercuryList);

  double writeParameter = (endtime/dt)/1000.0;

  sunMercuryClassical.calculateCenterofMass();
  sunMercuryClassical.writeForEach(writeParameter);
  sunMercuryClassical.writetoFile("sun_mercury/classical");

  sunMercuryRelativistic.calculateCenterofMass();
  sunMercuryRelativistic.relativistic();
  sunMercuryRelativistic.writeForEach(writeParameter);
  sunMercuryRelativistic.writetoFile("sun_mercury/relativistic");
  sunMercuryRelativistic.solveVelocityVerlet(endtime, dt);
}

/*
void solveMercuryPrecession(vector<Planet> &sunMercuryList, double endtime, double dt){

  System sunMercuryClassical("Sun-Mercury classical system", sunMercuryList);
  System sunMercuryRelativistic("Sun-Mercury relativistic system", sunMercuryList);

  sunMercuryClassical.calculateCenterofMass();
  sunMercuryClassical.writePerihelion("sun_mercury/classical");
  cout << 1 << endl;
  sunMercuryClassical.solveVelocityVerlet(endtime, dt);
  cout << 2 << endl;

  sunMercuryRelativistic.calculateCenterofMass();
  sunMercuryRelativistic.relativistic();
  sunMercuryRelativistic.writePerihelion("sun_mercury/relativistic");
  cout << 3 << endl;
  sunMercuryRelativistic.solveVelocityVerlet(endtime, dt);
  cout << 4 << endl;
}
*/
