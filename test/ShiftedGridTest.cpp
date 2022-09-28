#include "ShiftedGrid.h"
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
  // Construct a simple grid with shifted values

  auto shift = 2;
  NFmiDataMatrix<float> values(4, 2, kFloatMissing);
  for (auto j = 0UL; j < values.NY(); j++)
    for (auto i = 0UL; i < values.NX(); i++)
      values[i][j] = (i + shift) % values.NX();

  // 5x4 grid
  Fmi::CoordinateMatrix coords(5, 2, 0, 0, 1, 1);
  for (auto j = 0UL; j < coords.height(); j++)
    for (auto i = 0UL; i < coords.width(); i++)
      coords.set(i, j, (i - shift) % values.NX(), j);

  // Analyze validity
  auto analysis = std::make_shared<Fmi::CoordinateAnalysis>(Fmi::analysis(coords));

  // Contouring API
  SmartMet::Engine::Contour::ShiftedGrid grid(values, coords, analysis->valid, shift);

  Trax::Contour contour;
  contour.validate(true);

#if 1
  {
    Trax::IsobandLimits limits;
    limits.add(0, 1);

    auto results = contour.isobands(grid, limits);

    if (results.size() != 1)
      TEST_FAILED("Failed to contour valid area");

    string expected = "POLYGON ((0 0,0 1,1 1,1 0,0 0))";
    auto wkt = results[0].wkt();
    if (wkt != expected)
      TEST_FAILED("Incorrect result: " + wkt + " <> " + expected);
  }
  {
    Trax::IsobandLimits limits;
    limits.add(1, 2);

    auto results = contour.isobands(grid, limits);

    if (results.size() != 1)
      TEST_FAILED("Failed to contour valid area");

    string expected = "POLYGON ((1 0,1 1,2 1,2 0,1 0))";
    auto wkt = results[0].wkt();
    if (wkt != expected)
      TEST_FAILED("Incorrect result: " + wkt + " <> " + expected);
  }
#endif
#if 1
  {
    Trax::IsobandLimits limits;
    limits.add(0, 2);

    auto results = contour.isobands(grid, limits);

    if (results.size() != 1)
      TEST_FAILED("Failed to contour valid area");

    string expected = "POLYGON ((0 0,0 1,1 1,2 1,2 0,1 0,0 0))";
    auto wkt = results[0].wkt();
    if (wkt != expected)
      TEST_FAILED("Incorrect result: " + wkt + " <> " + expected);
  }
#endif
#if 1
  {
    Trax::IsobandLimits limits;
    limits.add(0, 3);

    auto results = contour.isobands(grid, limits);

    if (results.size() != 1)
      TEST_FAILED("Failed to contour valid area");

    string expected = "POLYGON ((0 0,0 1,1 1,2 1,3 1,3 0,2 0,1 0,0 0))";
    auto wkt = results[0].wkt();
    if (wkt != expected)
      TEST_FAILED("Incorrect result: " + wkt + " <> " + expected);
  }
#endif
#if 1
  {
    Trax::IsobandLimits limits;
    limits.add(1, 3);

    auto results = contour.isobands(grid, limits);

    if (results.size() != 1)
      TEST_FAILED("Failed to contour valid area");

    string expected = "POLYGON ((1 0,1 1,2 1,3 1,3 0,2 0,1 0))";
    auto wkt = results[0].wkt();
    if (wkt != expected)
      TEST_FAILED("Incorrect result: " + wkt + " <> " + expected);
  }
#endif
  TEST_PASSED();
}

class tests : public tframe::tests
{
  // Overridden message separator
  virtual const char *error_message_prefix() const { return "\n\t"; }
  // Main test suite
  void test() { TEST(normal); }
};
}  // namespace Tests

int main()
{
  cout << "\nShiftedGrid tester\n==================\n";
  Tests::tests t;
  return t.run();
}
