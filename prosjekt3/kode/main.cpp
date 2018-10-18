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

void taskC(Planet *, double, double);
void taskD(Planet *, double, double);
void taskE(Planet *, double, double);
void taskF(Planet *, double, double);
void taskG(Planet *, double, double);

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
  Planet *allplanets = new Planet[9];
  allplanets[0] = sun;
  for (int i = 0; i < 8; i++){
    allplanets[i+1] = extract(filename, i);
  }

  Planet *sunEarthList = new Planet[2];
  sunEarthList[0] = sun;
  sunEarthList[1] = earth;

  Planet *sunEarthJupiterList = new Planet[3];
  sunEarthJupiterList[0] = sun;
  sunEarthJupiterList[1] = earth;
  sunEarthJupiterList[2] = jupiter;

  Planet *sunMercuryList = new Planet[2];
  sunMercuryList[0] = sun;
  sunMercuryList[1] = mercury;
  if(scenario == 0){taskC(sunEarthList, endtime, dt);}
  if(scenario == 1){taskD(sunEarthList, endtime, dt);}
  if(scenario == 2){taskE(sunEarthJupiterList, endtime, dt);}
  if(scenario == 3){taskF(allplanets, endtime, dt);}
  if(scenario == 4){taskG(sunMercuryList, endtime, dt);}
  return 0;
}

void taskC(Planet *sunEarthList, double endtime, double dt){

  System sunEarth("Sun-Earth system", sunEarthList, 2);

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

void taskD(Planet *sunEarthList, double endtime, double dt){

  //flagg lese skaleringene fra kommandolinja
  double velosityScale = 2;
  double beta = 3.0;

  System sunEarthScale("Sun-Earth system", sunEarthList, 2);

  sunEarthScale.setBeta(beta);                               //endrer beta i kraftfunksjonen
  sunEarthScale.writetoFile("change_beta");
  sunEarthScale.solveVelocityVerlet(endtime, dt);
  sunEarthScale.setBeta(2);                                  //resetter beat til 2

  sunEarthScale.scalePlanetInitVel(velosityScale, 1);        //skalerer hastigheten til Jorda
  sunEarthScale.writetoFile("escape_velocity");
  sunEarthScale.solveVelocityVerlet(endtime, dt);
  double scaledVelocity = sunEarthScale.getPlanetInitVel(1); //henter skalert initsialhastighet

  cout << scaledVelocity << endl;
}

void taskE(Planet *sunEarthJupiterList, double endtime, double dt){

  System sunEarthJupiter("Sun-Earth_jupiter system", sunEarthJupiterList, 3);

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

void taskF(Planet *allplanets, double endtime, double dt){

  System solarsystem("Solar system", allplanets, 9);
  solarsystem.calculateCenterofMass();
  solarsystem.writetoFile("solarsystem");
  solarsystem.solveVelocityVerlet(endtime, dt);
}

void taskG(Planet *sunMercuryList, double endtime, double dt){

  System sunMercury("Sun-Mercury system", sunMercuryList, 2);
  sunMercury.writetoFile("sun_mercury/classical");
  sunMercury.solveVelocityVerlet(endtime, dt);

  sunMercury.relativistic();
  sunMercury.writetoFile("sun_mercury/relativistic");
  sunMercury.solveVelocityVerlet(endtime, dt);
}
