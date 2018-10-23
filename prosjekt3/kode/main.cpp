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
void varyVelocityBeta(vector<Planet>&, double, double, double, double);
void solveEarthJupiter(vector<Planet>&, double, double);
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


  Planet earth = extract(filename, 2);
  Planet jupiter = extract(filename, 4);
  Planet mercury = extract(filename, 0);

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
    double beta = atof(argv[5]);
    varyVelocityBeta(sunEarthList, endtime, dt, velocityscale, beta);
  }
  if(scenario == 3){solveEarthJupiter(sunEarthJupiterList, endtime, dt);}
  if(scenario == 4){solveSolarSystem(allplanets, endtime, dt);}
  if(scenario == 5){solveMercuryPrecession(sunMercuryList, endtime, dt);}

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
  sunEarth.solveForwardEuler(endtime, dt);

  sunEarth.writetoFile("euler_vs_verlet/verlet");
  sunEarth.solveVelocityVerlet(endtime, dt);
}

void varyVelocityBeta(vector<Planet>& sunEarthList, double endtime, double dt, double velocityScale, double beta){

  System sunEarthScale("Sun-Earth system", sunEarthList);

  sunEarthScale.setBeta(beta);                               //endrer beta i kraftfunksjonen
  sunEarthScale.writetoFile("change_beta/beta" + to_string(beta));
  sunEarthScale.solveVelocityVerlet(endtime, dt);
  sunEarthScale.setBeta(2);                                  //resetter beta til 2

  sunEarthScale.scalePlanetInitVel(velocityScale, 1);        //skalerer hastigheten til Jorda
  sunEarthScale.writetoFile("escape_velocity/velocityScale" + to_string(velocityScale));
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

void solveSolarSystem(vector<Planet> &allplanets, double endtime, double dt){

  System solarsystem("Solar system", allplanets);
  solarsystem.calculateCenterofMass();
  solarsystem.writetoFile("solarsystem");
  solarsystem.solveVelocityVerlet(endtime, dt);
}

void solveMercuryPrecession(vector<Planet> &sunMercuryList, double endtime, double dt){

  System sunMercury("Sun-Mercury system", sunMercuryList);
  sunMercury.writetoFile("sun_mercury/classical");
  sunMercury.solveVelocityVerlet(endtime, dt);

  sunMercury.relativistic();
  sunMercury.writetoFile("sun_mercury/relativistic");
  sunMercury.solveVelocityVerlet(endtime, dt);
}
