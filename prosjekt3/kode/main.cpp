
#include "coordinate.hpp"

using namespace std;

int main(){
  Coordinate a(2, 1, 2);
  Coordinate b(1, 2, 3);
  double r = a.norm();
  cout << r << endl;
  return 0;
}
