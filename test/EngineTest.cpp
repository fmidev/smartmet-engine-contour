#include "Engine.h"
#include <boost/date_time/posix_time/ptime.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/timer/timer.hpp>
#include <engines/querydata/Engine.h>
#include <gis/Box.h>
#include <gis/OGR.h>
#include <regression/tframe.h>
#include <spine/Options.h>
#include <spine/ParameterFactory.h>
#include <spine/Reactor.h>
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
  boost::posix_time::ptime t = boost::posix_time::time_from_string("2008-08-06 12:00");
  Spine::Parameter temperature = Spine::ParameterFactory::instance().parse("Temperature");
  q->param(temperature.number());

  // Full data area

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
    boost::hash_combine(valueshash, opt.data_hash_value());
    if (opt.level)
      q->selectLevel(*opt.level);

    auto matrix = qengine->getValues(q, valueshash, opt.time);
    auto geom = *(contour->contour(qhash, crs, crs, *matrix, *coords, opt).begin());

    auto result = Fmi::OGR::exportToSvg(*geom, Box::identity(), 1);
    string ok = "";
    if (result != ok)
      TEST_FAILED("Isovalue: 0\n\tExpected: " + ok + "\n\tObtained: " + result);
  }

  // above the maximum we get nothing
  {
    std::vector<double> isovalues{30};
    Engine::Contour::Options opt(temperature, t, isovalues);

    auto valueshash = qhash;
    boost::hash_combine(valueshash, opt.data_hash_value());
    if (opt.level)
      q->selectLevel(*opt.level);

    auto matrix = qengine->getValues(q, valueshash, opt.time);
    auto geom = *(contour->contour(qhash, crs, crs, *matrix, *coords, opt).begin());

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
    boost::hash_combine(valueshash, opt.data_hash_value());

    if (opt.level)
      q->selectLevel(*opt.level);

    auto matrix = qengine->getValues(q, valueshash, opt.time);
    auto geom = *(contour->contour(qhash, crs, crs, *matrix, *coords, opt).begin());

    auto result = Fmi::OGR::exportToSvg(*geom, Box::identity(), 1);

    string ok =
        "M5.7 52 5.9 52 6 51.9 6.2 51.9 6.3 51.9 6.4 51.8 6.5 51.8 6.7 51.7 6.9 51.6 7 51.4M15.7 "
        "52.3 15.9 52.3 16 52.3";
    if (result != ok)
      TEST_FAILED("Isovalue: 25\n\tExpected: " + ok + "\n\tObtained: " + result);

    // test another resolution
    result = Fmi::OGR::exportToSvg(*geom, Box::identity(), 2);

    ok = "M5.72 51.99 5.92 51.95 5.95 51.95 6.19 51.90 6.28 51.87 6.44 51.82 6.54 51.78 6.70 51.70 "
         "6.75 51.68 6.90 51.57 6.99 51.45M15.75 52.27 15.93 52.29 16 52.28";
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
    boost::hash_combine(valueshash, opt.data_hash_value());
    if (opt.level)
      q->selectLevel(*opt.level);

    auto matrix = qengine->getValues(q, valueshash, opt.time);
    auto geom = *(contour->contour(qhash, crs, crs, *matrix, *coords, opt).begin());

    auto result = Fmi::OGR::exportToSvg(*geom, Box::identity(), 1);
    string ok =
        "M5.7 52 5.9 52 6 51.9 6.2 51.9 6.3 51.9 6.4 51.8 6.5 51.8 6.7 51.7 6.9 51.6 7 51.5 7.1 "
        "51.5";
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
  boost::posix_time::ptime t = boost::posix_time::time_from_string("2008-08-06 12:00");
  Spine::Parameter temperature = Spine::ParameterFactory::instance().parse("Temperature");
  q->param(temperature.number());

  // Full data area

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
    boost::hash_combine(valueshash, opt.data_hash_value());
    if (opt.level)
      q->selectLevel(*opt.level);

    auto matrix = qengine->getValues(q, valueshash, opt.time);
    auto geom = *(contour->contour(qhash, crs, crs, *matrix, *coords, opt).begin());

    auto result = Fmi::OGR::exportToSvg(*geom, Box::identity(), 1);
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
    auto geom = *(contour->contour(qhash, crs, crs, *matrix, *coords, opt).begin());

    auto result = Fmi::OGR::exportToSvg(*geom, Box::identity(), 1);
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
    auto geom = *(contour->contour(qhash, crs, crs, *matrix, *coords, opt).begin());

    auto result = Fmi::OGR::exportToSvg(*geom, Box::identity(), 1);
    string ok =
        "M5.8 51.8 5.9 51.8 6 51.7 6.1 51.7 6.2 51.7 6.3 51.6 6.4 51.6 6.5 51.5 6.6 51.5 6.7 51.4 "
        "6.6 51.4 6.4 51.4 6.2 51.3 6 51.3 5.9 51.4 5.9 51.6 5.8 51.7Z";
    if (result != ok)
      TEST_FAILED("Expected: " + ok + "\n\tObtained: " + result);

    // test another resolution too
    result = Fmi::OGR::exportToSvg(*geom, Box::identity(), 2);
    ok = "M5.81 51.78 5.92 51.76 5.93 51.75 6.04 51.73 6.07 51.72 6.10 51.71 6.19 51.68 6.30 51.63 "
         "6.31 51.62 6.32 51.61 6.43 51.57 6.44 51.56 6.53 51.51 6.55 51.51 6.56 51.50 6.65 51.44 "
         "6.70 51.41 6.60 51.39 6.40 51.36 6.20 51.33 6 51.30 5.95 51.42 5.90 51.55 5.85 51.68Z";

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
    auto geom = *(contour->contour(qhash, crs, crs, *matrix, *coords, opt).begin());

    auto result = Fmi::OGR::exportToSvg(*geom, Box::identity(), 1);
    string ok =
        "M5.8 51.8 5.9 51.7 6 51.7 6.1 51.7 6.2 51.7 6.3 51.6 6.4 51.6 6.5 51.5 6.6 51.5 6.6 51.4 "
        "6.7 51.4 6.6 51.4 6.4 51.4 6.2 51.3 6 51.3 5.9 51.4 5.9 51.6 5.8 51.7Z";

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
      boost::hash_combine(valueshash, opt.data_hash_value());
      if (opt.level)
        q->selectLevel(*opt.level);

      auto matrix = qengine->getValues(q, valueshash, opt.time);
      auto geoms = contour->contour(qhash, crs, crs, *matrix, *coords, opt);
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
    boost::hash_combine(valueshash, opt.data_hash_value());
    if (opt.level)
      q->selectLevel(*opt.level);

    auto matrix = qengine->getValues(q, valueshash, opt.time);
    auto geoms = contour->contour(qhash, crs, crs, *matrix, *coords, opt);
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

  // Use native coordinates
  auto crs = q->SpatialReference();
  CoordinatesPtr coords = qengine->getWorldCoordinates(q);

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
      auto geoms = contour->contour(qhash, crs, crs, *matrix, *coords, opt);
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
    boost::hash_combine(valueshash, opt.data_hash_value());
    if (opt.level)
      q->selectLevel(*opt.level);

    auto matrix = qengine->getValues(q, valueshash, opt.time);
    auto geoms = contour->contour(qhash, crs, crs, *matrix, *coords, opt);
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

  // Use native coordinates
  auto crs = q->SpatialReference();
  CoordinatesPtr coords = qengine->getWorldCoordinates(q);

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
        auto geoms = contour->contour(qhash, crs, crs, *matrix, *coords, opt);
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

  // Use native coordinates
  auto crs = q->SpatialReference();
  CoordinatesPtr coords = qengine->getWorldCoordinates(q);

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
        auto geoms = contour->contour(qhash, crs, crs, *matrix, *coords, opt);
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

  Fmi::SpatialReference ykj("EPSG:2393");

  // Full data area

  auto world1 = q->area().XYToWorldXY(q->area().BottomLeft());
  auto world2 = q->area().XYToWorldXY(q->area().TopRight());
  Box area(world1.X(), world1.Y(), world2.X(), world2.Y(), 100, 100);
  std::size_t qhash = Engine::Querydata::hash_value(q);
  auto crs = q->SpatialReference();

  std::vector<Engine::Contour::Range> limits;
  limits.push_back(Engine::Contour::Range(-100, 0));
  Engine::Contour::Options opt(temperature, t, limits);  // freezing temperatures

  auto valueshash = qhash;
  boost::hash_combine(valueshash, opt.data_hash_value());
  if (opt.level)
    q->selectLevel(*opt.level);

  auto matrix = qengine->getValues(q, valueshash, opt.time);
  CoordinatesPtr coords = qengine->getWorldCoordinates(q, ykj);
  auto geom = *(contour->contour(qhash, crs, ykj, *matrix, *coords, opt).begin());

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

  // Use native coordinates
  auto crs = q->SpatialReference();
  CoordinatesPtr coords = qengine->getWorldCoordinates(q);

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

  auto geoms = contour->contour(qhash, crs, crs, *matrix, *coords, opt);

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
  reactor.init();
  qengine = reinterpret_cast<SmartMet::Engine::Querydata::Engine *>(
      reactor.getSingleton("Querydata", NULL));
  contour =
      reinterpret_cast<SmartMet::Engine::Contour::Engine *>(reactor.getSingleton("Contour", NULL));

  cout << endl << "Engine tester" << endl << "=============" << endl;
  Tests::tests t;
  return t.run();
}
