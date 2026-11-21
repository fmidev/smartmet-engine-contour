// Options::transform converts the contour limits given in converted units (value * multiplier +
// offset) back to the units of the data.

#include "Options.h"
#include <macgyver/DateTime.h>
#include <regression/tframe.h>
#include <spine/Parameter.h>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

namespace Tests
{
bool close(double value, double expected)
{
  return std::abs(value - expected) < 1e-9;
}

SmartMet::Spine::Parameter temperature()
{
  return {"Temperature", SmartMet::Spine::Parameter::Type::Data, kFmiTemperature};
}

// ----------------------------------------------------------------------

void transform_isovalues()
{
  using namespace SmartMet;
  Fmi::DateTime t = Fmi::DateTime::from_iso_string("20080805T030000");

  // Celsius to Fahrenheit: 50 F = 10 C, 68 F = 20 C
  Engine::Contour::Options both(temperature(), t, std::vector<double>{50, 68});
  both.transformation(1.8, 32);
  both.transform();
  if (!close(both.isovalues[0], 10) || !close(both.isovalues[1], 20))
    TEST_FAILED("Fahrenheit isovalues 50,68 should be 10,20 in Celsius, got " +
                std::to_string(both.isovalues[0]) + "," + std::to_string(both.isovalues[1]));

  // Kelvin to Celsius: offset only
  Engine::Contour::Options offset(temperature(), t, std::vector<double>{0});
  offset.offset = -273.15;
  offset.transform();
  if (!close(offset.isovalues[0], 273.15))
    TEST_FAILED("Celsius isovalue 0 should be 273.15 in Kelvin");

  // m/s to knots: multiplier only
  Engine::Contour::Options multiplier(temperature(), t, std::vector<double>{19.4384});
  multiplier.multiplier = 1.94384;
  multiplier.transform();
  if (!close(multiplier.isovalues[0], 10))
    TEST_FAILED("19.4384 knots should be 10 m/s");

  TEST_PASSED();
}

// ----------------------------------------------------------------------

void transform_limits()
{
  using namespace SmartMet;
  Fmi::DateTime t = Fmi::DateTime::from_iso_string("20080805T030000");

  std::vector<Engine::Contour::Range> limits;
  limits.emplace_back(50.0, 68.0);
  limits.emplace_back(std::optional<double>(68.0), std::optional<double>());

  Engine::Contour::Options opt(temperature(), t, limits);
  opt.transformation(1.8, 32);
  opt.transform();

  if (!close(*opt.limits[0].lolimit, 10) || !close(*opt.limits[0].hilimit, 20))
    TEST_FAILED("Fahrenheit limits 50...68 should be 10...20 in Celsius, got " +
                std::to_string(*opt.limits[0].lolimit) + "..." +
                std::to_string(*opt.limits[0].hilimit));
  if (!close(*opt.limits[1].lolimit, 20) || opt.limits[1].hilimit)
    TEST_FAILED("Open ended limit 68... should be 20... in Celsius");

  // A negative multiplier reverses the limits: -20...-10 becomes 10...20, and ...-10 becomes 10...
  std::vector<Engine::Contour::Range> neglimits;
  neglimits.emplace_back(-20.0, -10.0);
  neglimits.emplace_back(std::optional<double>(), std::optional<double>(-10.0));
  Engine::Contour::Options neg(temperature(), t, neglimits);
  neg.multiplier = -1;
  neg.transform();
  if (!close(*neg.limits[0].lolimit, 10) || !close(*neg.limits[0].hilimit, 20))
    TEST_FAILED("Limits -20...-10 with multiplier -1 should be 10...20");
  if (!neg.limits[1].lolimit || !close(*neg.limits[1].lolimit, 10) || neg.limits[1].hilimit)
    TEST_FAILED("Limit ...-10 with multiplier -1 should be 10...");

  TEST_PASSED();
}

// Test driver
class tests : public tframe::tests
{
  // Overridden message separator
  virtual const char* error_message_prefix() const { return "\n\t"; }
  // Main test suite
  void test()
  {
    TEST(transform_isovalues);
    TEST(transform_limits);
  }

};  // class tests

}  // namespace Tests

int main()
{
  cout << endl
       << "Engine::Contour::Options tester" << endl
       << "===============================" << endl;
  Tests::tests t;
  return t.run();
}
