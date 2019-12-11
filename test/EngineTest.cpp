#include "Engine.h"
#include <gis/Box.h>
#include <gis/OGR.h>
#include <regression/tframe.h>

#include <engines/querydata/Engine.h>
#include <spine/Options.h>
#include <spine/ParameterFactory.h>
#include <spine/Reactor.h>

#include <boost/date_time/posix_time/ptime.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/timer/timer.hpp>

#include <libconfig.h++>

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
  boost::posix_time::ptime t = boost::posix_time::time_from_string("2008-08-06 12:00");
  Spine::Parameter temperature = Spine::ParameterFactory::instance().parse("Temperature");
  q->param(temperature.number());

  // Full data area
  auto world1 = q->area().XYToWorldXY(q->area().BottomLeft());
  auto world2 = q->area().XYToWorldXY(q->area().TopRight());
  Box area(world1.X(), world1.Y(), world2.X(), world2.Y(), 100, 100);

  std::size_t qhash = Engine::Querydata::hash_value(q);
  std::string wkt = q->area().WKT();
  OGRSpatialReference *sr = nullptr;

  // Temperature for 200808061200 UTC:
  // Min:6.01 Mean:14.84 Max:25.95

  // below the minimum we get nothing
  {
    std::vector<double> isovalues;
    isovalues.push_back(0.0);
    Engine::Contour::Options opt(temperature, t, isovalues);

    auto valueshash = qhash;
    boost::hash_combine(valueshash, opt.data_hash_value());
    if (opt.level)
      q->selectLevel(*opt.level);

    auto matrix = qengine->getValues(q, valueshash, opt.time);
    CoordinatesPtr coords = qengine->getWorldCoordinates(q, sr);
    auto geom =
        *(contour->contour(qhash, wkt, *matrix, coords, opt, q->needsWraparound(), sr).begin());

    auto result = Fmi::OGR::exportToSvg(*geom, area, 1);
    string ok = "";
    if (result != ok)
      TEST_FAILED("Expected: " + ok + "\n\tObtained: " + result);
  }

  // above the maximum we get nothing
  {
    std::vector<double> isovalues;
    isovalues.push_back(30.0);
    Engine::Contour::Options opt(temperature, t, isovalues);

    auto valueshash = qhash;
    boost::hash_combine(valueshash, opt.data_hash_value());
    if (opt.level)
      q->selectLevel(*opt.level);

    auto matrix = qengine->getValues(q, valueshash, opt.time);
    CoordinatesPtr coords = qengine->getWorldCoordinates(q, sr);
    auto geom =
        *(contour->contour(qhash, wkt, *matrix, coords, opt, q->needsWraparound(), sr).begin());

    auto result = Fmi::OGR::exportToSvg(*geom, area, 1);
    string ok = "";
    if (result != ok)
      TEST_FAILED("Expected: " + ok + "\n\tObtained: " + result);
  }

  // should get something smallish just below the max
  {
    std::vector<double> isovalues;
    isovalues.push_back(25.5);
    Engine::Contour::Options opt(temperature, t, isovalues);

    auto valueshash = qhash;
    boost::hash_combine(valueshash, opt.data_hash_value());

    if (opt.level)
      q->selectLevel(*opt.level);

    auto matrix = qengine->getValues(q, valueshash, opt.time);
    CoordinatesPtr coords = qengine->getWorldCoordinates(q, sr);
    auto geom =
        *(contour->contour(qhash, wkt, *matrix, coords, opt, q->needsWraparound(), sr).begin());

    auto result = Fmi::OGR::exportToSvg(*geom, area, 1);
    string ok = "M0 97.4 0.7 97.9 0.9 98 1.5 98.6 1.6 98.7 2.2 99.3 2.2 99.4 2.6 100";
    if (result != ok)
      TEST_FAILED("Expected: " + ok + "\n\tObtained: " + result);

    // test another resolution
    result = Fmi::OGR::exportToSvg(*geom, area, 0);

    ok = "M0 97 1 98 1 99 2 99 3 100";
    if (result != ok)
      TEST_FAILED("Expected: " + ok + "\n\tObtained: " + result);
  }

  // Smoothen the data a little
  {
    std::vector<double> isovalues;
    isovalues.push_back(25.5);
    Engine::Contour::Options opt(temperature, t, isovalues);
    opt.filter_size = 1;
    opt.filter_degree = 1;

    auto valueshash = qhash;
    boost::hash_combine(valueshash, opt.data_hash_value());
    if (opt.level)
      q->selectLevel(*opt.level);

    auto matrix = qengine->getValues(q, valueshash, opt.time);
    CoordinatesPtr coords = qengine->getWorldCoordinates(q, sr);
    auto geom =
        *(contour->contour(qhash, wkt, *matrix, coords, opt, q->needsWraparound(), sr).begin());

    auto result = Fmi::OGR::exportToSvg(*geom, area, 1);
    string ok = "M0 97.5 0.7 97.9 0.8 98 1.5 98.6 1.6 98.7 2.1 99.3 2.2 99.5 2.6 100";
    if (result != ok)
      TEST_FAILED("Expected: " + ok + "\n\tObtained: " + result);
  }

  TEST_PASSED();
}

// ----------------------------------------------------------------------

void fills()
{
  using namespace SmartMet;
  using Fmi::Box;

  auto q = qengine->get("pal_skandinavia");
  boost::posix_time::ptime t = boost::posix_time::time_from_string("2008-08-06 12:00");
  Spine::Parameter temperature = Spine::ParameterFactory::instance().parse("Temperature");
  q->param(temperature.number());

  // Full data area

  auto world1 = q->area().XYToWorldXY(q->area().BottomLeft());
  auto world2 = q->area().XYToWorldXY(q->area().TopRight());
  Box area(world1.X(), world1.Y(), world2.X(), world2.Y(), 100, 100);
  std::size_t qhash = Engine::Querydata::hash_value(q);
  std::string wkt = q->area().WKT();
  OGRSpatialReference *sr = nullptr;

  // Temperature for 200808061200 UTC:
  // Min:6.01 Mean:14.84 Max:25.95

  // below the minimum we get nothing
  {
    std::vector<Engine::Contour::Range> limits;
    limits.push_back(Engine::Contour::Range(0.0, 5.0));
    Engine::Contour::Options opt(temperature, t, limits);

    auto valueshash = qhash;
    boost::hash_combine(valueshash, opt.data_hash_value());
    if (opt.level)
      q->selectLevel(*opt.level);

    auto matrix = qengine->getValues(q, valueshash, opt.time);
    CoordinatesPtr coords = qengine->getWorldCoordinates(q, sr);
    auto geom =
        *(contour->contour(qhash, wkt, *matrix, coords, opt, q->needsWraparound(), sr).begin());

    auto result = Fmi::OGR::exportToSvg(*geom, area, 1);
    string ok = "";
    if (result != ok)
      TEST_FAILED("Expected: " + ok + "\n\tObtained: " + result);
  }

  // above the maximum we get nothing
  {
    std::vector<Engine::Contour::Range> limits;
    limits.push_back(Engine::Contour::Range(30.0, 100.0));
    Engine::Contour::Options opt(temperature, t, limits);

    auto valueshash = qhash;
    boost::hash_combine(valueshash, opt.data_hash_value());
    if (opt.level)
      q->selectLevel(*opt.level);

    auto matrix = qengine->getValues(q, valueshash, opt.time);
    CoordinatesPtr coords = qengine->getWorldCoordinates(q, sr);
    auto geom =
        *(contour->contour(qhash, wkt, *matrix, coords, opt, q->needsWraparound(), sr).begin());

    auto result = Fmi::OGR::exportToSvg(*geom, area, 1);
    string ok = "";
    if (result != ok)
      TEST_FAILED("Expected: " + ok + "\n\tObtained: " + result);
  }

  // should get something smallish just below the max
  {
    std::vector<Engine::Contour::Range> limits;
    limits.push_back(Engine::Contour::Range(25.5, 30.0));
    Engine::Contour::Options opt(temperature, t, limits);

    auto valueshash = qhash;
    boost::hash_combine(valueshash, opt.data_hash_value());
    if (opt.level)
      q->selectLevel(*opt.level);

    auto matrix = qengine->getValues(q, valueshash, opt.time);
    CoordinatesPtr coords = qengine->getWorldCoordinates(q, sr);
    auto geom =
        *(contour->contour(qhash, wkt, *matrix, coords, opt, q->needsWraparound(), sr).begin());

    auto result = Fmi::OGR::exportToSvg(*geom, area, 1);
    string ok =
        "M0 100 0 99.3 0 98.7 0 98 0 97.4 0.3 97.6 0.4 97.6 0.7 97.9 0.8 97.9 0.9 98 1.2 98.3 1.5 "
        "98.6 1.6 98.7 1.9 99 2.2 99.3 2.2 99.4 2.5 99.8 2.6 100 2.2 100 1.5 100 0.7 100Z";
    if (result != ok)
      TEST_FAILED("Expected: " + ok + "\n\tObtained: " + result);

    // test another resolution too
    result = Fmi::OGR::exportToSvg(*geom, area, 0);
    ok = "M0 100 0 99 0 98 0 97 0 98 1 98 1 99 2 99 2 100 3 100 2 100 1 100Z";
    if (result != ok)
      TEST_FAILED("Expected: " + ok + "\n\tObtained: " + result);
  }

  // Smoothen the data a little
  {
    std::vector<Engine::Contour::Range> limits;
    limits.push_back(Engine::Contour::Range(25.5, 30.0));
    Engine::Contour::Options opt(temperature, t, limits);
    opt.filter_size = 1;
    opt.filter_degree = 1;

    auto valueshash = qhash;
    boost::hash_combine(valueshash, opt.data_hash_value());
    if (opt.level)
      q->selectLevel(*opt.level);

    auto matrix = qengine->getValues(q, valueshash, opt.time);
    CoordinatesPtr coords = qengine->getWorldCoordinates(q, sr);
    auto geom =
        *(contour->contour(qhash, wkt, *matrix, coords, opt, q->needsWraparound(), sr).begin());

    auto result = Fmi::OGR::exportToSvg(*geom, area, 1);
    string ok =
        "M0 100 0 99.3 0 98.7 0 98 0 97.5 0.3 97.7 0.5 97.8 0.7 97.9 0.8 98 1.2 98.3 1.5 98.6 1.6 "
        "98.7 1.9 99 2.1 99.3 2.2 99.4 2.2 99.5 2.5 99.8 2.6 100 2.2 100 1.5 100 0.7 100Z";
    if (result != ok)
      TEST_FAILED("Expected: " + ok + "\n\tObtained: " + result);
  }
  TEST_PASSED();
}

// ----------------------------------------------------------------------

void crossection()
{
  using namespace SmartMet;
  using Fmi::Box;

  auto q = qengine->get("hbm");
  boost::posix_time::ptime t = boost::posix_time::time_from_string("2014-07-28 02:00");
  Spine::Parameter temperature = Spine::ParameterFactory::instance().parse("TemperatureSea");
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
    string ok = "M54.5 25 59.8 30 59.8 25 74.7 22.5 59.8 20.3 55.5 24.3Z";
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
  boost::posix_time::ptime t = boost::posix_time::time_from_string("2015-03-13 12:00");
  Spine::Parameter temperature = Spine::ParameterFactory::instance().parse("Temperature");
  q->param(temperature.number());

  // Full data area

  auto world1 = q->area().XYToWorldXY(q->area().BottomLeft());
  auto world2 = q->area().XYToWorldXY(q->area().TopRight());
  Box area(world1.X(), world1.Y(), world2.X(), world2.Y(), 100, 100);
  std::size_t qhash = Engine::Querydata::hash_value(q);
  std::string wkt = q->area().WKT();
  OGRSpatialReference *sr = nullptr;

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
      boost::hash_combine(valueshash, opt.data_hash_value());
      if (opt.level)
        q->selectLevel(*opt.level);

      auto matrix = qengine->getValues(q, valueshash, opt.time);
      CoordinatesPtr coords = qengine->getWorldCoordinates(q, sr);
      auto geoms = contour->contour(qhash, wkt, *matrix, coords, opt, q->needsWraparound(), sr);
    }
  }
  TEST_PASSED();
}

void speed_all_at_once()
{
  using namespace SmartMet;
  using Fmi::Box;

  auto q = qengine->get("ecmwf_maailma_pinta");
  boost::posix_time::ptime t = boost::posix_time::time_from_string("2015-03-13 12:00");
  Spine::Parameter temperature = Spine::ParameterFactory::instance().parse("Temperature");
  q->param(temperature.number());

  // Full data area

  auto world1 = q->area().XYToWorldXY(q->area().BottomLeft());
  auto world2 = q->area().XYToWorldXY(q->area().TopRight());
  Box area(world1.X(), world1.Y(), world2.X(), world2.Y(), 100, 100);
  std::size_t qhash = Engine::Querydata::hash_value(q);
  std::string wkt = q->area().WKT();
  OGRSpatialReference *sr = nullptr;

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
    boost::hash_combine(valueshash, opt.data_hash_value());
    if (opt.level)
      q->selectLevel(*opt.level);

    auto matrix = qengine->getValues(q, valueshash, opt.time);
    CoordinatesPtr coords = qengine->getWorldCoordinates(q, sr);
    auto geoms = contour->contour(qhash, wkt, *matrix, coords, opt, q->needsWraparound(), sr);
  }
  TEST_PASSED();
}

// ----------------------------------------------------------------------

void pressure()
{
  using namespace SmartMet;
  using Fmi::Box;

  auto q = qengine->get("ecmwf_pressure");
  boost::posix_time::ptime t = boost::posix_time::time_from_string("2016-04-25 09:00");
  Spine::Parameter pressure = Spine::ParameterFactory::instance().parse("Pressure");
  q->param(pressure.number());

  // Full data area

  auto world1 = q->area().XYToWorldXY(q->area().BottomLeft());
  auto world2 = q->area().XYToWorldXY(q->area().TopRight());
  Box area(world1.X(), world1.Y(), world2.X(), world2.Y(), 360, 180);
  std::size_t qhash = Engine::Querydata::hash_value(q);
  std::string wkt = q->area().WKT();
  OGRSpatialReference *sr = nullptr;

  {
    std::cout << std::endl;
    boost::timer::auto_cpu_timer totaltimer(2, "\tAll contouring took %t sec CPU, %w sec real\n");

    for (int i = 950; i <= 1050; i += 5)
    {
      double value = i;

      std::string report = ("\tIsoline " + boost::lexical_cast<std::string>(value) +
                            " took %t sec CPU, %w sec real\n");
      boost::timer::auto_cpu_timer timer(2, report);

      std::vector<double> isolines;
      isolines.push_back(value);
      Engine::Contour::Options opt(pressure, t, isolines);

      auto valueshash = qhash;
      boost::hash_combine(valueshash, opt.data_hash_value());
      if (opt.level)
        q->selectLevel(*opt.level);

      auto matrix = qengine->getValues(q, valueshash, opt.time);
      CoordinatesPtr coords = qengine->getWorldCoordinates(q, sr);
      auto geoms = contour->contour(qhash, wkt, *matrix, coords, opt, q->needsWraparound(), sr);
    }
  }
  TEST_PASSED();
}

void pressure_all_at_once()
{
  using namespace SmartMet;
  using Fmi::Box;

  auto q = qengine->get("ecmwf_pressure");
  boost::posix_time::ptime t = boost::posix_time::time_from_string("2016-04-25 09:00");
  Spine::Parameter pressure = Spine::ParameterFactory::instance().parse("Pressure");
  q->param(pressure.number());

  // Full data area

  auto world1 = q->area().XYToWorldXY(q->area().BottomLeft());
  auto world2 = q->area().XYToWorldXY(q->area().TopRight());
  Box area(world1.X(), world1.Y(), world2.X(), world2.Y(), 360, 180);
  std::size_t qhash = Engine::Querydata::hash_value(q);
  std::string wkt = q->area().WKT();
  OGRSpatialReference *sr = nullptr;

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
    boost::hash_combine(valueshash, opt.data_hash_value());
    if (opt.level)
      q->selectLevel(*opt.level);

    auto matrix = qengine->getValues(q, valueshash, opt.time);
    CoordinatesPtr coords = qengine->getWorldCoordinates(q, sr);
    auto geoms = contour->contour(qhash, wkt, *matrix, coords, opt, q->needsWraparound(), sr);
  }
  TEST_PASSED();
}

// ----------------------------------------------------------------------

void fillvalidation()
{
  using namespace SmartMet;
  using Fmi::Box;

  auto q = qengine->get("ecmwf_temperature");
  Spine::Parameter temperature = Spine::ParameterFactory::instance().parse("Temperature");
  q->param(temperature.number());

  // Full data area

  auto world1 = q->area().XYToWorldXY(q->area().BottomLeft());
  auto world2 = q->area().XYToWorldXY(q->area().TopRight());
  Box area(world1.X(), world1.Y(), world2.X(), world2.Y(), 100, 100);
  std::size_t qhash = Engine::Querydata::hash_value(q);
  std::string wkt = q->area().WKT();
  OGRSpatialReference *sr = nullptr;

  {
    std::cout << std::endl;
    boost::timer::auto_cpu_timer totaltimer(2, "\tAll contouring took %t sec CPU, %w sec real\n");

    for (q->resetTime(); q->nextTime();)
    {
      // q->time(NFmiMetTime(2015,6,1,3,0,0));
      std::cout << "Time: " << q->validTime() << std::endl;
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
        Engine::Contour::Options opt(temperature, q->validTime(), limits);

        auto valueshash = qhash;
        boost::hash_combine(valueshash, opt.data_hash_value());
        if (opt.level)
          q->selectLevel(*opt.level);

        auto matrix = qengine->getValues(q, valueshash, opt.time);
        CoordinatesPtr coords = qengine->getWorldCoordinates(q, sr);
        auto geoms = contour->contour(qhash, wkt, *matrix, coords, opt, q->needsWraparound(), sr);
      }
    }
  }
  TEST_PASSED();
}

// ----------------------------------------------------------------------

void linevalidation()
{
  using namespace SmartMet;
  using Fmi::Box;

  auto q = qengine->get("ecmwf_temperature");
  Spine::Parameter temperature = Spine::ParameterFactory::instance().parse("Temperature");
  q->param(temperature.number());

  // Full data area

  auto world1 = q->area().XYToWorldXY(q->area().BottomLeft());
  auto world2 = q->area().XYToWorldXY(q->area().TopRight());
  Box area(world1.X(), world1.Y(), world2.X(), world2.Y(), 100, 100);
  std::size_t qhash = Engine::Querydata::hash_value(q);
  std::string wkt = q->area().WKT();
  OGRSpatialReference *sr = nullptr;

  {
    std::cout << std::endl;
    boost::timer::auto_cpu_timer totaltimer(2, "\tAll contouring took %t sec CPU, %w sec real\n");

    for (q->resetTime(); q->nextTime();)
    {
      // q->time(NFmiMetTime(2015,5,31,18,0,0));
      std::cout << "Time: " << q->validTime() << std::endl;
      for (int i = -50; i < 50; i += 2)
      // for(int i=-46; i<-44; i+=2)
      {
        double value = i;

        std::string report = ("\tIsoline " + boost::lexical_cast<std::string>(value) +
                              " took %t sec CPU, %w sec real\n");
        boost::timer::auto_cpu_timer timer(2, report);

        std::vector<double> isovalues;
        isovalues.push_back(value);
        Engine::Contour::Options opt(temperature, q->validTime(), isovalues);

        auto valueshash = qhash;
        boost::hash_combine(valueshash, opt.data_hash_value());
        if (opt.level)
          q->selectLevel(*opt.level);

        auto matrix = qengine->getValues(q, valueshash, opt.time);
        CoordinatesPtr coords = qengine->getWorldCoordinates(q, sr);
        auto geoms = contour->contour(qhash, wkt, *matrix, coords, opt, q->needsWraparound(), sr);
      }
    }
  }
  TEST_PASSED();
}

// ----------------------------------------------------------------------

void globalykj()
{
  using namespace SmartMet;
  using Fmi::Box;

  auto q = qengine->get("ecmwf_maailma_pinta");
  boost::posix_time::ptime t = boost::posix_time::time_from_string("2015-03-13 12:00");
  Spine::Parameter temperature = Spine::ParameterFactory::instance().parse("Temperature");
  q->param(temperature.number());

  // YKJ spatial reference

  std::unique_ptr<OGRSpatialReference> ykj(new OGRSpatialReference);
  ykj->SetFromUserInput("EPSG:2393");

  // Full data area

  auto world1 = q->area().XYToWorldXY(q->area().BottomLeft());
  auto world2 = q->area().XYToWorldXY(q->area().TopRight());
  Box area(world1.X(), world1.Y(), world2.X(), world2.Y(), 100, 100);
  std::size_t qhash = Engine::Querydata::hash_value(q);
  std::string wkt = q->area().WKT();
  OGRSpatialReference *sr = nullptr;

  std::vector<Engine::Contour::Range> limits;
  limits.push_back(Engine::Contour::Range(-100, 0));
  Engine::Contour::Options opt(temperature, t, limits);  // freezing temperatures

  auto valueshash = qhash;
  boost::hash_combine(valueshash, opt.data_hash_value());
  if (opt.level)
    q->selectLevel(*opt.level);

  auto matrix = qengine->getValues(q, valueshash, opt.time);
  CoordinatesPtr coords = qengine->getWorldCoordinates(q, sr);
  auto geom =
      *(contour->contour(qhash, wkt, *matrix, coords, opt, q->needsWraparound(), sr).begin());

  if (!geom)
    TEST_FAILED("Failed to contour temperature interval -100...0");

  OGREnvelope env;
  geom->getEnvelope(&env);

#if 0
	std::cout << std::fixed << std::setprecision(1)
			  << "\tX: " << env.MinX << "..." << env.MaxX
			  << "\tY: " << env.MinY << "..." << env.MaxY
			  << std::endl;
#endif

  if (env.MinX > -2670969)
    TEST_FAILED("Engine::Contour MinX should be < -2670969");
  if (env.MaxX < 11621694)
    TEST_FAILED("Engine::Contour MaxX should be > 11621694");
  if (env.MinY > -10003576)
    TEST_FAILED("Engine::Contour MinY should be < -10003576");
  if (env.MaxY < 10445033)
    TEST_FAILED("Engine::Contour MaxY should be > 10445033");
  TEST_PASSED();
}

// ----------------------------------------------------------------------

void worldwrap()
{
  using namespace SmartMet;
  using Fmi::Box;

  auto q = qengine->get("gfs");
  q->firstTime();
  boost::posix_time::ptime t = q->validTime();
  Spine::Parameter temperature = Spine::ParameterFactory::instance().parse("Temperature");
  q->param(temperature.number());

  // Full data area

  auto world1 = q->area().XYToWorldXY(q->area().BottomLeft());
  auto world2 = q->area().XYToWorldXY(q->area().TopRight());
  Box area(world1.X(), world1.Y(), world2.X(), world2.Y(), 100, 100);
  std::size_t qhash = Engine::Querydata::hash_value(q);
  std::string wkt = q->area().WKT();

  OGRSpatialReference *sr = nullptr;

  // This contour spans the world horizontally
  double lolimit = 0;
  double hilimit = 2;

  std::vector<Engine::Contour::Range> limits;
  limits.push_back(Engine::Contour::Range(lolimit, hilimit));
  Engine::Contour::Options opt(temperature, t, limits);

  auto valueshash = qhash;
  boost::hash_combine(valueshash, opt.data_hash_value());
  if (opt.level)
    q->selectLevel(*opt.level);

  auto matrix = qengine->getValues(q, valueshash, opt.time);
  CoordinatesPtr coords = qengine->getWorldCoordinates(q, sr);

  auto geoms = contour->contour(qhash, wkt, *matrix, coords, opt, q->needsWraparound(), sr);

  if (geoms.empty())
    TEST_FAILED("Failed to contour GFS data interval 0-2");

  OGREnvelope envelope;
  geoms[0]->getEnvelope(&envelope);
  if (envelope.MinX != 0)
    TEST_FAILED("Contour 0-2 minimum x value should be 0, not " + std::to_string(envelope.MinX));
  if (envelope.MaxX != 360)
    TEST_FAILED("Contour 0-2 maximum x value should be 360, not " + std::to_string(envelope.MaxX));

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
#if 1
    TEST(lines);
    contour->clearCache();
    TEST(fills);
    contour->clearCache();
    TEST(crossection);
    contour->clearCache();
    TEST(speed);
    contour->clearCache();
    TEST(speed_all_at_once);
    contour->clearCache();
    TEST(pressure);
    contour->clearCache();
    TEST(pressure_all_at_once);
    contour->clearCache();
    TEST(worldwrap);
    contour->clearCache();
#else
    // these have been used only to make sure everything validates. Too slow for other testing
    TEST(fillvalidation);
    TEST(linevalidation);
    TEST(globalykj);
#endif
  }

};  // class tests

}  // namespace Tests

int main(void)
{
  SmartMet::Spine::Options opts;
  opts.configfile = "cnf/reactor.conf";
  opts.parseConfig();

  SmartMet::Spine::Reactor reactor(opts);
  qengine = reinterpret_cast<SmartMet::Engine::Querydata::Engine *>(
      reactor.getSingleton("Querydata", NULL));
  contour =
      reinterpret_cast<SmartMet::Engine::Contour::Engine *>(reactor.getSingleton("Contour", NULL));

  cout << endl << "Engine tester" << endl << "=============" << endl;
  Tests::tests t;
  t.run();
}
