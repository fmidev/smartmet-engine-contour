#include "Engine.h"
#include <boost/lexical_cast.hpp>
#include <boost/timer/timer.hpp>
#include <engines/querydata/Engine.h>
#include <gis/Box.h>
#include <gis/OGR.h>
#include <macgyver/DateTime.h>
#include <macgyver/Hash.h>
#include <macgyver/StringConversion.h>
#include <regression/tframe.h>
#include <spine/Options.h>
#include <spine/Reactor.h>
#include <timeseries/ParameterFactory.h>
#include <libconfig.h++>
#include <ogr_api.h>
#include <ogr_geometry.h>
#include <ogr_spatialref.h>

using namespace std;

std::shared_ptr<SmartMet::Engine::Querydata::Engine> qengine;
std::shared_ptr<SmartMet::Engine::Contour::Engine> contour;

// Note: pal_skandinavia = stereographic,20,90,60:6,51.3,49,70.2

namespace Tests
{
// ----------------------------------------------------------------------

void lines()
{
  using namespace SmartMet;
  using Fmi::Box;

  auto q = qengine->get("pal_skandinavia");
  Fmi::DateTime t = Fmi::DateTime::from_string("2008-08-06 12:00");
  Spine::Parameter temperature = TimeSeries::ParameterFactory::instance().parse("Temperature");
  q->param(temperature.number());

  // Full data area
  auto world1 = q->area().XYToWorldXY(q->area().BottomLeft());
  auto world2 = q->area().XYToWorldXY(q->area().TopRight());
  Box area(world1.X(), world1.Y(), world2.X(), world2.Y(), 100, 100);

  std::size_t qhash = Engine::Querydata::hash_value(q);

  // Use native coordinates
  auto crs = q->SpatialReference();
  CoordinatesPtr coords = qengine->getWorldCoordinates(q);

  // Temperature for 200808061200 UTC:
  // Min:6.01 Mean:14.84 Max:25.95

  // below the minimum we get nothing
  {
    std::vector<double> isovalues{0};
    Engine::Contour::Options opt(temperature, t, isovalues);

    auto valueshash = qhash;
    Fmi::hash_combine(valueshash, opt.data_hash_value());
    if (opt.level)
      q->selectLevel(*opt.level);

    auto matrix = qengine->getValues(q, valueshash, opt.time);
    auto geom = *(contour->contour(qhash, crs, *matrix, *coords, opt).begin());

    auto result = Fmi::OGR::exportToSvg(*geom, area, 1);
    string ok = "";
    if (result != ok)
      TEST_FAILED("Isovalue: 0\n\tExpected: " + ok + "\n\tObtained: " + result);
  }

  // above the maximum we get nothing
  {
    std::vector<double> isovalues{30};
    Engine::Contour::Options opt(temperature, t, isovalues);

    auto valueshash = qhash;
    Fmi::hash_combine(valueshash, opt.data_hash_value());
    if (opt.level)
      q->selectLevel(*opt.level);

    auto matrix = qengine->getValues(q, valueshash, opt.time);
    auto geom = *(contour->contour(qhash, crs, *matrix, *coords, opt).begin());

    auto result = Fmi::OGR::exportToSvg(*geom, Box::identity(), 1);
    string ok = "";
    if (result != ok)
      TEST_FAILED("Isovalue: 30\n\tExpected: " + ok + "\n\tObtained: " + result);
  }

  // should get something smallish just below the max
  {
    std::vector<double> isovalues{25};
    Engine::Contour::Options opt(temperature, t, isovalues);

    auto valueshash = qhash;
    Fmi::hash_combine(valueshash, opt.data_hash_value());

    if (opt.level)
      q->selectLevel(*opt.level);

    auto matrix = qengine->getValues(q, valueshash, opt.time);
    auto geom = *(contour->contour(qhash, crs, *matrix, *coords, opt).begin());

    auto result = Fmi::OGR::exportToSvg(*geom, area, 1);

    string ok =
        "M0 96.3 0.6 96.6 0.7 96.7 1.5 97.1 1.8 97.3 2.2 97.7 2.5 98 3 98.5 3.1 98.7 3.5 99.3 3.7 "
        "100M35.2 100 35.8 99.9 36.1 100";

    if (result != ok)
      TEST_FAILED("Isovalue: 25\n\tExpected: " + ok + "\n\tObtained: " + result);

    // test another resolution
    result = Fmi::OGR::exportToSvg(*geom, area, 2);

    ok = "M0 96.31 0.63 96.64 0.75 96.71 1.49 97.12 1.78 97.32 2.24 97.7 2.55 97.99 2.99 98.5 3.12 "
         "98.66 3.5 99.33 3.68 100M35.17 100 35.82 99.93 36.06 100";
    if (result != ok)
      TEST_FAILED("Isovalue: 25\n\tExpected: " + ok + "\n\tObtained: " + result);
  }

  // Smoothen the data a little
  {
    std::vector<double> isovalues{25};
    Engine::Contour::Options opt(temperature, t, isovalues);
    opt.filter_size = 1;
    opt.filter_degree = 1;

    auto valueshash = qhash;
    Fmi::hash_combine(valueshash, opt.data_hash_value());
    if (opt.level)
      q->selectLevel(*opt.level);

    auto matrix = qengine->getValues(q, valueshash, opt.time);
    auto geom = *(contour->contour(qhash, crs, *matrix, *coords, opt).begin());

    auto result = Fmi::OGR::exportToSvg(*geom, area, 1);
    string ok =
        "M0 96.3 0.7 96.6 0.7 96.7 1.5 97.1 1.7 97.3 2.2 97.7 2.5 98 3 98.5 3.1 98.7 3.6 99.3 3.7 "
        "99.6 3.9 100";
    if (result != ok)
      TEST_FAILED("Isovalue: 25 smoothed\n\tExpected: " + ok + "\n\tObtained: " + result);
  }

  TEST_PASSED();
}

// ----------------------------------------------------------------------

void fills()
{
  using namespace SmartMet;
  using Fmi::Box;

  auto q = qengine->get("pal_skandinavia");
  Fmi::DateTime t = Fmi::DateTime::from_string("2008-08-06 12:00");
  Spine::Parameter temperature = TimeSeries::ParameterFactory::instance().parse("Temperature");
  q->param(temperature.number());

  // Full data area
  auto world1 = q->area().XYToWorldXY(q->area().BottomLeft());
  auto world2 = q->area().XYToWorldXY(q->area().TopRight());
  Box area(world1.X(), world1.Y(), world2.X(), world2.Y(), 100, 100);
  std::size_t qhash = Engine::Querydata::hash_value(q);

  // Use native coordinates
  auto crs = q->SpatialReference();
  CoordinatesPtr coords = qengine->getWorldCoordinates(q);

  // Temperature for 200808061200 UTC:
  // Min:6.01 Mean:14.84 Max:25.95

  // below the minimum we get nothing
  {
    std::vector<Engine::Contour::Range> limits;
    limits.push_back(Engine::Contour::Range(0.0, 5.0));
    Engine::Contour::Options opt(temperature, t, limits);

    auto valueshash = qhash;
    Fmi::hash_combine(valueshash, opt.data_hash_value());
    if (opt.level)
      q->selectLevel(*opt.level);

    auto matrix = qengine->getValues(q, valueshash, opt.time);
    auto geom = *(contour->contour(qhash, crs, *matrix, *coords, opt).begin());

    auto result = Fmi::OGR::exportToSvg(*geom, area, 1);
    string ok = "";
    if (result != ok)
      TEST_FAILED("Isoband: 0-5\n\tExpected: " + ok + "\n\tObtained: " + result);
  }

  // above the maximum we get nothing
  {
    std::vector<Engine::Contour::Range> limits;
    limits.push_back(Engine::Contour::Range(30.0, 100.0));
    Engine::Contour::Options opt(temperature, t, limits);

    auto valueshash = qhash;
    Fmi::hash_combine(valueshash, opt.data_hash_value());
    if (opt.level)
      q->selectLevel(*opt.level);

    auto matrix = qengine->getValues(q, valueshash, opt.time);
    auto geom = *(contour->contour(qhash, crs, *matrix, *coords, opt).begin());

    auto result = Fmi::OGR::exportToSvg(*geom, area, 1);
    string ok = "";
    if (result != ok)
      TEST_FAILED("Isoband: 30-100\n\tExpected: " + ok + "\n\tObtained: " + result);
  }

  // should get something smallish just below the max
  {
    std::vector<Engine::Contour::Range> limits;
    limits.push_back(Engine::Contour::Range(25.5, 30.0));
    Engine::Contour::Options opt(temperature, t, limits);

    auto valueshash = qhash;
    Fmi::hash_combine(valueshash, opt.data_hash_value());
    if (opt.level)
      q->selectLevel(*opt.level);

    auto matrix = qengine->getValues(q, valueshash, opt.time);
    auto geom = *(contour->contour(qhash, crs, *matrix, *coords, opt).begin());

    auto result = Fmi::OGR::exportToSvg(*geom, area, 1);
    string ok =
        "M0 100 0 99.3 0 98.7 0 98 0 97.4 0.7 97.9 0.9 98 1.5 98.6 1.6 98.7 2.2 99.3 2.2 99.4 2.6 "
        "100 2.2 100 1.5 100 0.7 100Z";

    if (result != ok)
      TEST_FAILED("Isoband (1 decimal): 25.5-30\n\tExpected: " + ok + "\n\tObtained: " + result);

    // test another resolution too
    result = Fmi::OGR::exportToSvg(*geom, area, 2);
    ok = "M0 100 0 99.33 0 98.66 0 97.99 0 97.42 0.75 97.87 0.92 97.99 1.49 98.58 1.56 98.66 2.16 "
         "99.33 2.24 99.42 2.6 100 2.24 100 1.49 100 0.75 100Z";
    if (result != ok)
      TEST_FAILED("Isoband (2 decimals): 25.5-30 \n\tExpected: " + ok + "\n\tObtained: " + result);
  }

  // Smoothen the data a little
  {
    std::vector<Engine::Contour::Range> limits;
    limits.push_back(Engine::Contour::Range(25.5, 30.0));
    Engine::Contour::Options opt(temperature, t, limits);
    opt.filter_size = 1;
    opt.filter_degree = 1;

    auto valueshash = qhash;
    Fmi::hash_combine(valueshash, opt.data_hash_value());
    if (opt.level)
      q->selectLevel(*opt.level);

    auto matrix = qengine->getValues(q, valueshash, opt.time);
    auto geom = *(contour->contour(qhash, crs, *matrix, *coords, opt).begin());

    auto result = Fmi::OGR::exportToSvg(*geom, area, 1);
    string ok =
        "M0 100 0 99.3 0 98.7 0 98 0 97.4 0.7 97.9 0.8 98 1.5 98.6 1.6 98.7 2.1 99.3 2.2 99.4 2.6 "
        "100 2.2 100 1.5 100 0.7 100Z";

    if (result != ok)
      TEST_FAILED("Smoothened isoband: 25.5-30\n\tExpected: " + ok + "\n\tObtained: " + result);
  }
  TEST_PASSED();
}

// ----------------------------------------------------------------------

// Exercise the Trax grid smoother path (Options::smoother), independent of the
// legacy Savitzky-Golay filter_size/filter_degree path.

void trax_smoother()
{
  using namespace SmartMet;
  using Fmi::Box;

  auto q = qengine->get("pal_skandinavia");
  Fmi::DateTime t = Fmi::DateTime::from_string("2008-08-06 12:00");
  Spine::Parameter temperature = TimeSeries::ParameterFactory::instance().parse("Temperature");
  q->param(temperature.number());

  auto world1 = q->area().XYToWorldXY(q->area().BottomLeft());
  auto world2 = q->area().XYToWorldXY(q->area().TopRight());
  Box area(world1.X(), world1.Y(), world2.X(), world2.Y(), 100, 100);
  std::size_t qhash = Engine::Querydata::hash_value(q);

  auto crs = q->SpatialReference();
  CoordinatesPtr coords = qengine->getWorldCoordinates(q);

  auto svg_isoband = [&](const std::optional<Trax::SmoothOptions> &smoother) -> string
  {
    std::vector<Engine::Contour::Range> limits;
    limits.push_back(Engine::Contour::Range(25.5, 30.0));
    Engine::Contour::Options opt(temperature, t, limits);
    opt.smoother = smoother;

    auto valueshash = qhash;
    Fmi::hash_combine(valueshash, opt.data_hash_value());
    if (opt.level)
      q->selectLevel(*opt.level);

    auto matrix = qengine->getValues(q, valueshash, opt.time);
    auto geom = *(contour->contour(qhash, crs, *matrix, *coords, opt).begin());
    return Fmi::OGR::exportToSvg(*geom, area, 1);
  };

  // A box smoother (radius 1, single pass) must change the contour, and the
  // result is anchored as a regression golden value.
  {
    Trax::SmoothOptions box;
    box.method = Trax::SmoothMethod::Box;
    box.radius = 1;
    box.passes = 1;

    auto plain = svg_isoband(std::nullopt);
    auto smoothed = svg_isoband(box);

    if (smoothed == plain)
      TEST_FAILED("Box smoother did not change the isoband 25.5-30:\n\t" + smoothed);

    string ok =
        "M0 100 0 99.3 0 98.7 0 98 0 97.7 0.7 97.9 0.8 98 1.5 98.6 1.6 98.7 2.1 99.3 2.2 99.6 2.4 "
        "100 2.2 100 1.5 100 0.7 100Z";
    if (smoothed != ok)
      TEST_FAILED("Box-smoothed isoband 25.5-30\n\tExpected: " + ok + "\n\tObtained: " + smoothed);
  }

  // An inactive smoother (radius 0) must be a no-op: identical to no smoother.
  {
    Trax::SmoothOptions inactive;
    inactive.method = Trax::SmoothMethod::Box;
    inactive.radius = 0;
    if (svg_isoband(inactive) != svg_isoband(std::nullopt))
      TEST_FAILED("Inactive box smoother changed the result");
  }

  // The median smoother must also run and produce a non-empty contour.
  {
    Trax::SmoothOptions median;
    median.method = Trax::SmoothMethod::Median;
    median.radius = 1;
    median.passes = 1;
    if (svg_isoband(median).empty())
      TEST_FAILED("Median smoother produced an empty isoband 25.5-30");
  }

  TEST_PASSED();
}

// ----------------------------------------------------------------------

void crossection()
{
  using namespace SmartMet;
  using Fmi::Box;

  auto q = qengine->get("hbm");
  Fmi::DateTime t = Fmi::DateTime::from_string("2014-07-28 02:00");
  Spine::Parameter temperature = TimeSeries::ParameterFactory::instance().parse("TemperatureSea");
  q->param(temperature.number());

  // Full data area

  Box area = Box::identity();

  double lon1 = 24.9;
  double lat1 = 60.2;  // Helsinki
  double lon2 = 24.7;
  double lat2 = 59.4;  // Tallinna
  std::size_t steps = 3;

  // Min:1.80 Mean:10.15 Max:24.60
  // At the bottom of the Finnish Gulf it is about 4-5 degrees,
  // at the surface around 20 degrees at best

  // with a large isoband we get the sea cross section
  {
    std::vector<Engine::Contour::Range> limits;
    limits.push_back(Engine::Contour::Range(0, 30));
    Engine::Contour::Options opt(temperature, t, limits);

    std::shared_ptr<NFmiFastQueryInfo> qInfo = q->info();

    auto geom = *(contour->crossection(*qInfo, opt, lon1, lat1, lon2, lat2, steps).begin());

    auto result = Fmi::OGR::exportToSvg(*geom, area, 1);
    string ok =
        "M29.9 0 29.9 5 29.9 10 29.9 15 29.9 20 29.9 25 59.8 30 59.8 25 89.7 20 89.7 15 89.7 10 "
        "89.7 5 89.7 0 59.8 0Z";

    if (result != ok)
      TEST_FAILED("Expected: " + ok + "\n\tObtained: " + result);
  }

  // just above minimum should be at the bottom of the sea
  {
    std::vector<Engine::Contour::Range> limits;
    limits.push_back(Engine::Contour::Range(0, 5.0));
    Engine::Contour::Options opt(temperature, t, limits);
    std::shared_ptr<NFmiFastQueryInfo> qInfo = q->info();

    auto geom = *(contour->crossection(*qInfo, opt, lon1, lat1, lon2, lat2, steps).begin());

    auto result = Fmi::OGR::exportToSvg(*geom, area, 1);
    string ok = "M54.5 25 59.8 30 59.8 25 74.7 22.5 59.8 20.3Z";
    if (result != ok)
      TEST_FAILED("Expected: " + ok + "\n\tObtained: " + result);
  }
  TEST_PASSED();
}

// ----------------------------------------------------------------------

void speed()
{
  using namespace SmartMet;
  using Fmi::Box;

  auto q = qengine->get("ecmwf_maailma_pinta");
  Fmi::DateTime t = Fmi::DateTime::from_string("2015-03-13 12:00");
  Spine::Parameter temperature = TimeSeries::ParameterFactory::instance().parse("Temperature");
  q->param(temperature.number());

  // Full data area

  auto world1 = q->area().XYToWorldXY(q->area().BottomLeft());
  auto world2 = q->area().XYToWorldXY(q->area().TopRight());
  Box area(world1.X(), world1.Y(), world2.X(), world2.Y(), 100, 100);
  std::size_t qhash = Engine::Querydata::hash_value(q);

  // Use native coordinates
  auto crs = q->SpatialReference();
  CoordinatesPtr coords = qengine->getWorldCoordinates(q);

  {
    std::cout << std::endl;
    boost::timer::auto_cpu_timer totaltimer(2, "\tAll contouring took %t sec CPU, %w sec real\n");

    for (int i = -50; i < 50; i += 2)
    {
      double lolimit = i;
      double hilimit = i + 2;

      std::string report = ("\tIsoband " + Fmi::to_string(lolimit) + "..." +
                            Fmi::to_string(hilimit) + " took %t sec CPU, %w sec real\n");
      boost::timer::auto_cpu_timer timer(2, report);

      std::vector<Engine::Contour::Range> limits;
      limits.push_back(Engine::Contour::Range(lolimit, hilimit));
      Engine::Contour::Options opt(temperature, t, limits);

      auto valueshash = qhash;
      Fmi::hash_combine(valueshash, opt.data_hash_value());
      if (opt.level)
        q->selectLevel(*opt.level);

      auto matrix = qengine->getValues(q, valueshash, opt.time);
      auto geoms = contour->contour(qhash, crs, *matrix, *coords, opt);
    }
  }
  TEST_PASSED();
}

void speed_all_at_once()
{
  using namespace SmartMet;
  using Fmi::Box;

  auto q = qengine->get("ecmwf_maailma_pinta");
  Fmi::DateTime t = Fmi::DateTime::from_string("2015-03-13 12:00");
  Spine::Parameter temperature = TimeSeries::ParameterFactory::instance().parse("Temperature");
  q->param(temperature.number());

  // Full data area

  auto world1 = q->area().XYToWorldXY(q->area().BottomLeft());
  auto world2 = q->area().XYToWorldXY(q->area().TopRight());
  Box area(world1.X(), world1.Y(), world2.X(), world2.Y(), 100, 100);
  std::size_t qhash = Engine::Querydata::hash_value(q);

  // Use native coordinates
  auto crs = q->SpatialReference();
  CoordinatesPtr coords = qengine->getWorldCoordinates(q);

  {
    std::cout << std::endl;
    boost::timer::auto_cpu_timer totaltimer(2, "\tAll contouring took %t sec CPU, %w sec real\n");

    std::vector<Engine::Contour::Range> limits;
    for (int i = -50; i < 50; i += 2)
    {
      double lolimit = i;
      double hilimit = i + 2;
      limits.push_back(Engine::Contour::Range(lolimit, hilimit));
    }

    Engine::Contour::Options opt(temperature, t, limits);

    auto valueshash = qhash;
    Fmi::hash_combine(valueshash, opt.data_hash_value());
    if (opt.level)
      q->selectLevel(*opt.level);

    auto matrix = qengine->getValues(q, valueshash, opt.time);
    auto geoms = contour->contour(qhash, crs, *matrix, *coords, opt);
  }
  TEST_PASSED();
}

// ----------------------------------------------------------------------

void pressure()
{
  using namespace SmartMet;
  using Fmi::Box;

  auto q = qengine->get("ecmwf_pressure");
  Fmi::DateTime t = Fmi::DateTime::from_string("2016-04-25 09:00");
  Spine::Parameter pressure = TimeSeries::ParameterFactory::instance().parse("Pressure");
  q->param(pressure.number());

  // Full data area

  auto world1 = q->area().XYToWorldXY(q->area().BottomLeft());
  auto world2 = q->area().XYToWorldXY(q->area().TopRight());
  Box area(world1.X(), world1.Y(), world2.X(), world2.Y(), 360, 180);
  std::size_t qhash = Engine::Querydata::hash_value(q);

  // Use native coordinates
  auto crs = q->SpatialReference();
  CoordinatesPtr coords = qengine->getWorldCoordinates(q);

  {
    std::cout << std::endl;
    boost::timer::auto_cpu_timer totaltimer(2, "\tAll contouring took %t sec CPU, %w sec real\n");

    for (int i = 950; i <= 1050; i += 5)
    // for (int i = 1000; i <= 1000; i += 5)
    {
      double value = i;

      std::string report =
          ("\tIsoline " + Fmi::to_string(value) + " took %t sec CPU, %w sec real\n");
      boost::timer::auto_cpu_timer timer(2, report);

      std::vector<double> isolines;
      isolines.push_back(value);
      Engine::Contour::Options opt(pressure, t, isolines);

      auto valueshash = qhash;
      Fmi::hash_combine(valueshash, opt.data_hash_value());
      if (opt.level)
        q->selectLevel(*opt.level);

      auto matrix = qengine->getValues(q, valueshash, opt.time);
      auto geoms = contour->contour(qhash, crs, *matrix, *coords, opt);
    }
  }
  TEST_PASSED();
}

void pressure_all_at_once()
{
  using namespace SmartMet;
  using Fmi::Box;

  auto q = qengine->get("ecmwf_pressure");
  Fmi::DateTime t = Fmi::DateTime::from_string("2016-04-25 09:00");
  Spine::Parameter pressure = TimeSeries::ParameterFactory::instance().parse("Pressure");
  q->param(pressure.number());

  // Full data area

  auto world1 = q->area().XYToWorldXY(q->area().BottomLeft());
  auto world2 = q->area().XYToWorldXY(q->area().TopRight());
  Box area(world1.X(), world1.Y(), world2.X(), world2.Y(), 360, 180);
  std::size_t qhash = Engine::Querydata::hash_value(q);

  // Use native coordinates
  auto crs = q->SpatialReference();
  CoordinatesPtr coords = qengine->getWorldCoordinates(q);

  {
    std::cout << std::endl;
    boost::timer::auto_cpu_timer totaltimer(2, "\tAll contouring took %t sec CPU, %w sec real\n");

    std::vector<double> isolines;
    for (int i = 950; i <= 1050; i += 5)
    {
      isolines.push_back(i);
    }

    Engine::Contour::Options opt(pressure, t, isolines);

    auto valueshash = qhash;
    Fmi::hash_combine(valueshash, opt.data_hash_value());
    if (opt.level)
      q->selectLevel(*opt.level);

    auto matrix = qengine->getValues(q, valueshash, opt.time);
    auto geoms = contour->contour(qhash, crs, *matrix, *coords, opt);
  }
  TEST_PASSED();
}

// ----------------------------------------------------------------------

// Area of a possibly empty or missing geometry

// A tile to be contoured and the hash of the expected SVG representation of the
// result. The expected values have been captured from the implementation which
// tested every cell of the grid for an overlap with the clipping box without
// caching the mask of the contoured cells.

struct TileTest
{
  double x1;
  double y1;
  double x2;
  double y2;
  std::size_t expected;
};

// Contour the given tiles and compare the results with the expected ones. Also
// verify that caching the mask of the contoured cells does not alter the result:
// the mask is cached for each grid and clipping box, and is hence shared by
// requests for different isobands of the same tile.

void compare_tiled_contour(const std::string &theName,
                           std::size_t theHash,
                           const Fmi::SpatialReference &theCRS,
                           const NFmiDataMatrix<float> &theMatrix,
                           const Fmi::CoordinateMatrix &theCoordinates,
                           const SmartMet::Engine::Contour::Options &theOptions,
                           const SmartMet::Engine::Contour::Options &theOtherOptions,
                           const std::vector<TileTest> &theTiles)
{
  using namespace SmartMet;
  using Fmi::Box;

  for (const auto &t : theTiles)
  {
    contour->clearCache();

    Box tile(t.x1, t.y1, t.x2, t.y2, 100, 100);

    auto geom =
        contour->contour(theHash, theCRS, theMatrix, theCoordinates, tile, theOptions).at(0);

    const auto svg = Fmi::OGR::exportToSvg(*geom, tile, 1);
    const auto result = Fmi::hash_value(svg);

    if (result != t.expected)
      TEST_FAILED(theName + ": tile " + std::to_string(t.x1) + "," + std::to_string(t.y1) + "," +
                  std::to_string(t.x2) + "," + std::to_string(t.y2) +
                  "\n\tExpected hash: " + std::to_string(t.expected) +
                  "\n\tObtained hash: " + std::to_string(result) + " for " + svg);

    // The very same tile is contoured again for another isoband, which finds the
    // mask of the contoured cells from the cache. The result must not change.

    auto cold =
        contour->contour(theHash, theCRS, theMatrix, theCoordinates, tile, theOtherOptions).at(0);
    const auto cold_svg = Fmi::OGR::exportToSvg(*cold, tile, 1);

    contour->clearCache();

    auto warm =
        contour->contour(theHash, theCRS, theMatrix, theCoordinates, tile, theOtherOptions).at(0);
    const auto warm_svg = Fmi::OGR::exportToSvg(*warm, tile, 1);

    if (cold_svg != warm_svg)
      TEST_FAILED(theName + ": a cached cell mask changed the result\n\tExpected: " + warm_svg +
                  "\n\tObtained: " + cold_svg);
  }

  // Every cell of the grid overlaps a clipping box this large, and hence the
  // result must be identical to contouring without a clipping box at all.

  const double large = 1E10;
  Box everything(-large, -large, large, large, 100, 100);

  contour->clearCache();
  auto untiled = contour->contour(theHash, theCRS, theMatrix, theCoordinates, theOptions).at(0);
  const auto untiled_svg = Fmi::OGR::exportToSvg(*untiled, everything, 5);

  contour->clearCache();
  auto clipped =
      contour->contour(theHash, theCRS, theMatrix, theCoordinates, everything, theOptions).at(0);
  const auto clipped_svg = Fmi::OGR::exportToSvg(*clipped, everything, 5);

  if (untiled_svg != clipped_svg)
    TEST_FAILED(theName +
                ": contouring with a clipping box covering everything differs from contouring "
                "without a clipping box\n\tExpected: " +
                untiled_svg + "\n\tObtained: " + clipped_svg);
}

// Tiles covering the quadrants of the given box, plus a tile completely outside
// of it so that no cell of the grid overlaps the clipping box.

std::vector<TileTest> quadrant_tiles(
    double x1, double y1, double x2, double y2, const std::vector<std::size_t> &theExpected)
{
  std::vector<TileTest> ret;

  const auto xm = 0.5 * (x1 + x2);
  const auto ym = 0.5 * (y1 + y2);
  const auto w = x2 - x1;

  ret.push_back({x1, y1, xm, ym, 0});
  ret.push_back({xm, y1, x2, ym, 0});
  ret.push_back({x1, ym, xm, y2, 0});
  ret.push_back({xm, ym, x2, y2, 0});
  ret.push_back({x1 + 10 * w, y1, x1 + 11 * w, ym, 0});  // outside the data

  for (auto i = 0UL; i < ret.size() && i < theExpected.size(); i++)
    ret[i].expected = theExpected[i];

  return ret;
}

void tiles()
{
  using namespace SmartMet;

  // 1. A regular grid contoured in its native projection

  {
    auto q = qengine->get("pal_skandinavia");
    Fmi::DateTime t = Fmi::DateTime::from_string("2008-08-06 12:00");
    Spine::Parameter temperature = TimeSeries::ParameterFactory::instance().parse("Temperature");
    q->param(temperature.number());

    std::size_t qhash = Engine::Querydata::hash_value(q);
    auto crs = q->SpatialReference();
    CoordinatesPtr coords = qengine->getWorldCoordinates(q);

    std::vector<Engine::Contour::Range> limits{Engine::Contour::Range(10.0, 15.0)};
    Engine::Contour::Options opt(temperature, t, limits);

    std::vector<Engine::Contour::Range> other_limits{Engine::Contour::Range(15.0, 20.0)};
    Engine::Contour::Options other_opt(temperature, t, other_limits);

    auto valueshash = qhash;
    Fmi::hash_combine(valueshash, opt.data_hash_value());
    auto matrix = qengine->getValues(q, valueshash, opt.time);

    auto world1 = q->area().XYToWorldXY(q->area().BottomLeft());
    auto world2 = q->area().XYToWorldXY(q->area().TopRight());

    const std::vector<std::size_t> expected{1926903029129084511UL,
                                            1987795063613586975UL,
                                            17777160359691528594UL,
                                            1260149825643717534UL,
                                            6142509188972423790UL};  // empty result

    auto tiles = quadrant_tiles(world1.X(), world1.Y(), world2.X(), world2.Y(), expected);

    compare_tiled_contour("pal_skandinavia", qhash, crs, *matrix, *coords, opt, other_opt, tiles);
  }

  // 2. Global data reprojected to a stereographic projection. Coordinates
  //    outside the projection are missing, and the cells of the grid are both
  //    huge and invalid near the projection boundaries.

  {
    auto q = qengine->get("gfs");
    q->firstTime();
    Fmi::DateTime t = q->validTime();
    Spine::Parameter temperature = TimeSeries::ParameterFactory::instance().parse("Temperature");
    q->param(temperature.number());

    std::size_t qhash = Engine::Querydata::hash_value(q);

    auto pal = qengine->get("pal_skandinavia");
    auto crs = pal->SpatialReference();
    CoordinatesPtr coords = qengine->getWorldCoordinates(q, crs);

    std::vector<Engine::Contour::Range> limits{Engine::Contour::Range(0.0, 5.0)};
    Engine::Contour::Options opt(temperature, t, limits);

    std::vector<Engine::Contour::Range> other_limits{Engine::Contour::Range(5.0, 10.0)};
    Engine::Contour::Options other_opt(temperature, t, other_limits);

    auto valueshash = qhash;
    Fmi::hash_combine(valueshash, opt.data_hash_value());
    auto matrix = qengine->getValues(q, valueshash, opt.time);

    auto world1 = pal->area().XYToWorldXY(pal->area().BottomLeft());
    auto world2 = pal->area().XYToWorldXY(pal->area().TopRight());

    const std::vector<std::size_t> expected{15090222048163346963UL,
                                            9836952456016350948UL,
                                            15821242027694908769UL,
                                            5165353953160476314UL,
                                            6142509188972423790UL};  // empty result

    auto tiles = quadrant_tiles(world1.X(), world1.Y(), world2.X(), world2.Y(), expected);

    compare_tiled_contour(
        "gfs in stereographic", qhash, crs, *matrix, *coords, opt, other_opt, tiles);
  }

  TEST_PASSED();
}

// ----------------------------------------------------------------------

// Contouring many isobands of a single tile of a reprojected global grid. Each
// isoband misses the contour cache, but all of them share the mask of the cells
// to be contoured, which is calculated only once for the tile.

void tile_speed()
{
  using namespace SmartMet;
  using Fmi::Box;

  auto q = qengine->get("gfs");
  q->firstTime();
  Fmi::DateTime t = q->validTime();
  Spine::Parameter temperature = TimeSeries::ParameterFactory::instance().parse("Temperature");
  q->param(temperature.number());

  std::size_t qhash = Engine::Querydata::hash_value(q);
  auto pal = qengine->get("pal_skandinavia");
  auto crs = pal->SpatialReference();
  CoordinatesPtr coords = qengine->getWorldCoordinates(q, crs);

  auto world1 = pal->area().XYToWorldXY(pal->area().BottomLeft());
  auto world2 = pal->area().XYToWorldXY(pal->area().TopRight());

  // A small tile in the middle of the requested area
  const auto xm = 0.5 * (world1.X() + world2.X());
  const auto ym = 0.5 * (world1.Y() + world2.Y());
  const auto dx = 0.1 * (world2.X() - world1.X());
  const auto dy = 0.1 * (world2.Y() - world1.Y());
  Box tile(xm - dx, ym - dy, xm + dx, ym + dy, 256, 256);

  std::vector<Engine::Contour::Range> dummy{Engine::Contour::Range(0.0, 1.0)};
  Engine::Contour::Options dummy_opt(temperature, t, dummy);
  auto valueshash = qhash;
  Fmi::hash_combine(valueshash, dummy_opt.data_hash_value());
  auto matrix = qengine->getValues(q, valueshash, dummy_opt.time);

  std::cout << std::endl;
  {
    boost::timer::auto_cpu_timer timer(3, "\t50 isobands of one tile: %t sec CPU, %w sec real\n");

    for (int i = 0; i < 50; i++)
    {
      // Distinct limits so that the contour cache always misses
      std::vector<Engine::Contour::Range> limits{
          Engine::Contour::Range(-50.0 + 2 * i, -49.0 + 2 * i)};
      Engine::Contour::Options opt(temperature, t, limits);
      contour->contour(qhash, crs, *matrix, *coords, tile, opt);
    }
  }

  TEST_PASSED();
}

// ----------------------------------------------------------------------

void worldwrap()
{
  using namespace SmartMet;
  using Fmi::Box;

  auto q = qengine->get("gfs");
  q->firstTime();
  Fmi::DateTime t = q->validTime();
  Spine::Parameter temperature = TimeSeries::ParameterFactory::instance().parse("Temperature");
  q->param(temperature.number());

  // Full data area

  std::size_t qhash = Engine::Querydata::hash_value(q);

  // Use native coordinates
  auto crs = q->SpatialReference();
  CoordinatesPtr coords = qengine->getWorldCoordinates(q, "WGS84");

  // This contour spans the world horizontally
  double lolimit = 0;
  double hilimit = 2;

  std::vector<Engine::Contour::Range> limits;
  limits.push_back(Engine::Contour::Range(lolimit, hilimit));
  Engine::Contour::Options opt(temperature, t, limits);

  auto valueshash = qhash;
  Fmi::hash_combine(valueshash, opt.data_hash_value());
  if (opt.level)
    q->selectLevel(*opt.level);

  auto matrix = qengine->getValues(q, valueshash, opt.time);

  auto geoms = contour->contour(qhash, crs, *matrix, *coords, opt);

  if (geoms.empty())
    TEST_FAILED("Failed to contour GFS data interval 0-2");

  // GFS data is from 0 to 359.75, but getWorldCoordinates should return an extended matrix to 360.

  OGREnvelope envelope;
  geoms[0]->getEnvelope(&envelope);
  if (std::abs(envelope.MaxX - 180) > 0.01)
    TEST_FAILED("Contour 0-2 maximum x value should be 180, not " + std::to_string(envelope.MaxX));
  // if (std::abs(envelope.MinX - (-180)) > 0.01)
  if (std::abs(envelope.MinX - (-179.75)) > 0.01)
    TEST_FAILED("Contour 0-2 minimum x value should be -180, not " + std::to_string(envelope.MinX));

  TEST_PASSED();
}

// Test driver
class tests : public tframe::tests
{
  // Overridden message separator
  virtual const char *error_message_prefix() const { return "\n\t"; }
  // Main test suite
  void test()
  {
    TEST(lines);
    contour->clearCache();
    TEST(fills);
    contour->clearCache();
    TEST(trax_smoother);
    contour->clearCache();
    TEST(crossection);
    contour->clearCache();
    TEST(worldwrap);
    contour->clearCache();
    TEST(tiles);
    contour->clearCache();
    TEST(tile_speed);
    contour->clearCache();
    TEST(pressure);
    contour->clearCache();
    TEST(pressure_all_at_once);
    contour->clearCache();
    TEST(speed);
    contour->clearCache();
    TEST(speed_all_at_once);
    contour->clearCache();
  }

};  // class tests

}  // namespace Tests

int main(void)
{
  SmartMet::Spine::Options opts;
  opts.configfile = "cnf/reactor.conf";
  opts.parseConfig();

  SmartMet::Spine::Reactor reactor(opts);
  reactor.init();
  qengine = reactor.getEngine<SmartMet::Engine::Querydata::Engine>("Querydata", nullptr);
  contour = reactor.getEngine<SmartMet::Engine::Contour::Engine>("Contour", nullptr);

  cout << endl << "Engine tester" << endl << "=============" << endl;
  Tests::tests t;
  auto result = t.run();
  qengine.reset();
  contour.reset();
  reactor.shutdown();
  return result;
}
