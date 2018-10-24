#include "catch.hpp"
#include "coordinate.hpp"


TEST_CASE("COORDINATE: VECTOR ARITHMETIC"){
  Coordinate u(2, -1, 2);
  Coordinate v(3, 4, -2);

  Coordinate sumUV = u + v;
  Coordinate subUV = u - v;
  Coordinate lhsmulU = 2.0*u;
  Coordinate rhsmulU= u*2.0;
  Coordinate divU = u/3.0;

  Coordinate true_sumUV(5, 3, 0);
  Coordinate true_subUV(-1, -5, 4);
  Coordinate true_mulU(4, -2, 4);
  Coordinate true_divU(2.0/3.0, -1.0/3.0, 2.0/3.0);

  REQUIRE(sumUV == true_sumUV);
  REQUIRE(subUV == true_subUV);
  REQUIRE(lhsmulU == true_mulU);
  REQUIRE(rhsmulU == true_mulU);
  REQUIRE(divU == true_divU);
}

TEST_CASE("COORDINATE: NORM"){
  Coordinate u(2, -1, 2);

  double normU = u.norm();

  REQUIRE(normU == 3);
}

TEST_CASE("COORDINATE: CROSS PRODUCT"){
  Coordinate u(3, -5, 2);
  Coordinate v(-2, 1, 6);

  Coordinate uCrossv = u^v;

  Coordinate true_uCrossv(-32, -22, -7);

  REQUIRE(uCrossv== true_uCrossv);
}

TEST_CASE("COORDINATE: DOT PRODUCT"){
  Coordinate u(3, -5, 2);
  Coordinate v(-2, 1, 6);

  double uDotv = u*v;

  double true_uDotv = 1;

  REQUIRE(uDotv == uDotv);
}
