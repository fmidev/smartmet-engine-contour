#include "Engine.h"
#include <macgyver/DateTime.h>
#include <boost/lexical_cast.hpp>
#include <boost/timer/timer.hpp>
#include <engines/querydata/Engine.h>
#include <gis/Box.h>
#include <gis/OGR.h>
#include <macgyver/Hash.h>
#include <regression/tframe.h>
#include <spine/Options.h>
#include <spine/Reactor.h>
#include <timeseries/ParameterFactory.h>
#include <libconfig.h++>
#include <ogr_geometry.h>
#include <ogr_spatialref.h>

using namespace std;

SmartMet::Engine::Querydata::Engine *qengine;
SmartMet::Engine::Contour::Engine *contour;

// Note: pal_skandinavia = stereographic,20,90,60:6,51.3,49,70.2

namespace Tests
{
// ----------------------------------------------------------------------

void lines()
{
  using namespace SmartMet;
  using Fmi::Box;

  auto q = qengine->get("pal_skandinavia");
  Fmi::DateTime t = boost::posix_time::time_from_string("2008-08-06 12:00");
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
  Fmi::DateTime t = boost::posix_time::time_from_string("2008-08-06 12:00");
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

void crossection()
{
  using namespace SmartMet;
  using Fmi::Box;

  auto q = qengine->get("hbm");
  Fmi::DateTime t = boost::posix_time::time_from_string("2014-07-28 02:00");
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

    boost::shared_ptr<NFmiFastQueryInfo> qInfo = q->info();

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
    boost::shared_ptr<NFmiFastQueryInfo> qInfo = q->info();

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
  Fmi::DateTime t = boost::posix_time::time_from_string("2015-03-13 12:00");
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

      std::string report =
          ("\tIsoband " + boost::lexical_cast<std::string>(lolimit) + "..." +
           boost::lexical_cast<std::string>(hilimit) + " took %t sec CPU, %w sec real\n");
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
  Fmi::DateTime t = boost::posix_time::time_from_string("2015-03-13 12:00");
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
  Fmi::DateTime t = boost::posix_time::time_from_string("2016-04-25 09:00");
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

      std::string report = ("\tIsoline " + boost::lexical_cast<std::string>(value) +
                            " took %t sec CPU, %w sec real\n");
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
  Fmi::DateTime t = boost::posix_time::time_from_string("2016-04-25 09:00");
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
    TEST(crossection);
    contour->clearCache();
    TEST(worldwrap);
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
  qengine = reinterpret_cast<SmartMet::Engine::Querydata::Engine *>(
      reactor.getSingleton("Querydata", NULL));
  contour =
      reinterpret_cast<SmartMet::Engine::Contour::Engine *>(reactor.getSingleton("Contour", NULL));

  cout << endl << "Engine tester" << endl << "=============" << endl;
  Tests::tests t;
  return t.run();
}
