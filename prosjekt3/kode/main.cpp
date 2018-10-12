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

void taskC(Planet *);
void taskD(Planet *);
void taskE(Planet *);
void taskF(Planet *);
void taskG(Planet *);

int main(){
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

<<<<<<< HEAD
  taskC(sunEarthList);
  taskD(sunEarthList);
=======
  //taskC(sunEarthList);
  //taskD(sunEarthList);
>>>>>>> 62892cbfc12a8f0e9cadce7dcdd6ab22bcb5270d
  taskE(sunEarthJupiterList);
  taskF(allplanets);
  taskG(sunMercuryList);

}

void taskC(Planet *sunEarthList){

  //FLAGG lese inn endtime og dt fra kommandolinja
  double dt = 0.001;
  double endtime = 10;

  System sunEarth("Sun-Earth system", sunEarthList, 2);

  clock_t c_start = clock();
  sunEarth.solveEuler(endtime, dt);
  clock_t c_end = clock();
  sunEarth.writetoFile("../data/euler_vs_verlet/euler");

  // Beregner CPU-tid i milisekunder
  double eulerTime = 1000.0 * (c_end - c_start) / CLOCKS_PER_SEC;

  c_start = clock();
  sunEarth.solveVerlet(endtime, dt);
  c_end = clock();
  sunEarth.writetoFile("../data/euler_vs_verlet/verlet");

  // Beregner CPU-tid i milisekunder
  double verletTime = 1000.0 * (c_end - c_start) / CLOCKS_PER_SEC;

  cout << eulerTime << " " << verletTime << endl;
}

void taskD(Planet *sunEarthList){

  //flagg lese skaleringene fra kommandolinja
  double velosityScale = 2;
  double beta = 2.5;
  double dt = 0.001;
  double endtime = 2;

  System sunEarthScale("Sun-Earth system", sunEarthList, 2);

  sunEarthScale.setBeta(beta);                               //endrer beta i kraftfunksjonen
  sunEarthScale.solveVerlet(endtime, dt);
  sunEarthScale.writetoFile("../data/change_beta");
  sunEarthScale.setBeta(2);                                  //resetter beat til 2

  sunEarthScale.scalePlanetInitVel(velosityScale, 1);        //skalerer hastigheten til Jorda
  sunEarthScale.solveVerlet(endtime, dt);
  sunEarthScale.writetoFile("../data/escape_velocity");
  double scaledVelocity = sunEarthScale.getPlanetInitVel(1); //henter skalert initsialhastighet

  cout << scaledVelocity << endl;
}

void taskE(Planet *sunEarthJupiterList){
  //FLAGG lese inn dt og endtime fra kommandolinja
  double dt = 0.001;
  double endtime = 10;

  System sunEarthJupiter("Sun-Earth_jupiter system", sunEarthJupiterList, 3);

  //løser for Jupiters masse skalert med 1
  sunEarthJupiter.solveVerlet(endtime, dt);
  sunEarthJupiter.writetoFile("../data/sun_earth_jupiter/jupiter_mass_1");

  //løser for jupiters masse skalert med 10
  sunEarthJupiter.scalePlanetMass(10,2);
  sunEarthJupiter.solveVerlet(endtime, dt);
  sunEarthJupiter.writetoFile("../data/sun_earth_jupiter/jupiter_mass_10");

  //løser for jupiters masse skalert med 1000 = 10*100
  sunEarthJupiter.scalePlanetMass(100,2);
  sunEarthJupiter.solveVerlet(endtime, dt);
  sunEarthJupiter.writetoFile("../data/sun_earth_jupiter/jupiter_mass_1000");
}

void taskF(Planet *allplanets){
  //FLAGG lese inn dt og endtime fra kommandolinja
  double dt = 0.001;
  double endtime = 10;

  System solarsystem("Solar system", allplanets, 9);
  solarsystem.calculateCenterofMass();
  solarsystem.solveVerlet(endtime,dt);
  solarsystem.writetoFile("../data/solarsystem");
}

void taskG(Planet *sunMercuryList){
  //Flagg lese inn dt fra kommandolinja?
  double dt = 0.001;
  double endtime = 100;

  System sunMercury("Sun-Mercury system", sunMercuryList, 2);
  sunMercury.solveVerlet(endtime,dt);
  sunMercury.writetoFile("../data/sun_mercury/classical");

  sunMercury.relativistic("on");
  sunMercury.solveVerlet(endtime,dt);
  sunMercury.writetoFile("../data/sun_mercury/relativistic");
}
