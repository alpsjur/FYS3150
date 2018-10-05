
#include "coordinate.hpp"

using namespace std;

int main(){
  Coordinate pos(1, 1, 1);
  Coordinate vel(1, 1, 1);
  Coordinate acc(1, 1, 1);
  pos = pos + vel + acc;
  cout << pos;
  return 0;
}
