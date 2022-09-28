#include "Grid.h"
#include <gis/CoordinateMatrix.h>
#include <gis/CoordinateMatrixAnalysis.h>
#include <newbase/NFmiDataMatrix.h>
#include <regression/tframe.h>
#include <trax/Contour.h>
#include <trax/IsobandLimits.h>
#include <limits>

using namespace std;

namespace Tests
{
void normal()
{
  // Construct a simple grid with one corner value missing, and try contouring NaN...NaN
  // to see if we get the expected missing area

  NFmiDataMatrix<float> values(2, 2, kFloatMissing);
  values[0][0] = 0;
  values[0][1] = 2;
  values[1][0] = 4;
  values[1][1] = 6;

  // 2x2 grid from 0,0 to 10,10
  Fmi::CoordinateMatrix coords(2, 2, 0, 0, 10, 10);

  // Analyze validity
  auto analysis = std::make_shared<Fmi::CoordinateAnalysis>(Fmi::analysis(coords));

  // Contouring API
  SmartMet::Engine::Contour::Grid grid(values, coords, analysis->valid);
  grid.shell(50);  // external shell for missing values

  // Isoband for missing values
  Trax::IsobandLimits limits;
  limits.add(2, 4);

  Trax::Contour contour;
  contour.validate(true);
  auto results = contour.isobands(grid, limits);

  if (results.size() != 1)
    TEST_FAILED("Failed to contour valid area");

  string expected = "POLYGON ((0 10,5 10,10 0,5 0,0 10))";
  auto wkt = results[0].wkt();
  if (wkt != expected)
    TEST_FAILED("Incorrect result: " + wkt + " <> " + expected);

  TEST_PASSED();
}

void valid()
{
  // Construct a simple grid with one corner value missing, and try contouring NaN...NaN
  // to see if we get the expected missing area

  NFmiDataMatrix<float> values(2, 2, kFloatMissing);
  values[0][0] = 0;
  values[0][1] = 1;
  values[1][0] = 2;
  values[1][1] = 5;

  // 2x2 grid from 0,0 to 1,1
  Fmi::CoordinateMatrix coords(2, 2, 0, 0, 1, 1);

  // Analyze validity
  auto analysis = std::make_shared<Fmi::CoordinateAnalysis>(Fmi::analysis(coords));

  // Contouring API
  SmartMet::Engine::Contour::Grid grid(values, coords, analysis->valid);
  grid.shell(10);  // external shell for missing values

  // Isoband for missing values
  Trax::IsobandLimits limits;
  limits.add(-std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity());

  Trax::Contour contour;
  contour.validate(true);
  auto results = contour.isobands(grid, limits);

  if (results.size() != 1)
    TEST_FAILED("Failed to contour valid area");

  string expected = "POLYGON ((0 0,0 1,1 1,1 0,0 0))";
  auto wkt = results[0].wkt();
  if (wkt != expected)
    TEST_FAILED("Incorrect result: " + wkt + " <> " + expected);

  TEST_PASSED();
}

void missing()
{
  // Construct a simple grid with one corner value missing, and try contouring NaN...NaN
  // to see if we get the expected missing area

  NFmiDataMatrix<float> values(2, 2, kFloatMissing);
  values[0][0] = std::numeric_limits<float>::quiet_NaN();
  values[0][1] = 1;
  values[1][0] = 2;
  values[1][1] = 5;

  // 2x2 grid from -1,-1 to 1,1
  Fmi::CoordinateMatrix coords(2, 2, -1, -1, 1, 1);

  // Analyze validity
  auto analysis = std::make_shared<Fmi::CoordinateAnalysis>(Fmi::analysis(coords));

  // Contouring API
  SmartMet::Engine::Contour::Grid grid(values, coords, analysis->valid);
  grid.shell(2);  // external shell for missing values

  // Isoband for missing values
  Trax::IsobandLimits limits;
  limits.add(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN());

  Trax::Contour contour;
  contour.validate(true);
  auto results = contour.isobands(grid, limits);

  if (results.size() != 1)
    TEST_FAILED("Failed to contour missing area");

  string ok =
      "POLYGON ((-2 -2,-2 -1,-2 1,-2 2,-1 2,1 2,2 2,2 1,2 -1,2 -2,1 -2,-1 -2,-2 -2),(1 -1,1 1,-1 "
      "1,1 -1))";
  auto wkt = results[0].wkt();
  if (wkt != ok)
    TEST_FAILED("Incorrect result: " + wkt + " <> " + ok);

  TEST_PASSED();
}

class tests : public tframe::tests
{
  // Overridden message separator
  virtual const char *error_message_prefix() const { return "\n\t"; }
  // Main test suite
  void test()
  {
    TEST(normal);
    TEST(valid);
    TEST(missing);
  }
};
}  // namespace Tests

int main()
{
  cout << "\nGrid tester\n===========\n";
  Tests::tests t;
  return t.run();
}
