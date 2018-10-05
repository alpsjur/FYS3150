
#include "coordinate.hpp"

using namespace std;

int main(){
  Coordinate pos(2, 1, 2);
  Coordinate vel(1, 2, 3);
  double dt = 0.001;
  pos = pos + dt*vel;
  cout << pos;
  return 0;
}
