// ======================================================================
/*!
 * \brief Implementation details for the Contour engine
 */
// ======================================================================

#include "Impl.h"
#include "Options.h"
#include <boost/functional/hash.hpp>
#include <gdal/ogr_core.h>
#include <gdal/ogr_geometry.h>
#include <gdal/ogr_spatialref.h>
#include <geos/geom/GeometryFactory.h>
#include <geos/io/WKBWriter.h>
#include <macgyver/StringConversion.h>
#include <macgyver/WorkQueue.h>
#include <spine/Exception.h>
#include <tron/SavitzkyGolay2D.h>
#include <cmath>

// ----------------------------------------------------------------------
/*!
 * \brief Hash for OGR spatial reference
 *
 * Compiler lookup issues prevented placing this into spine.
 * Placing this into the above anonymous namespace or into the
 * Contour namespace below did not work either.
 */
// ----------------------------------------------------------------------

std::size_t hash_value(const OGRSpatialReference &theSR)
{
  try
  {
    char *wkt;
    theSR.exportToWkt(&wkt);
    std::string tmp(wkt);
    OGRFree(wkt);
    boost::hash<std::string> hasher;
    return hasher(tmp);
  }
  catch (...)
  {
    throw SmartMet::Spine::Exception(BCP, "Operation failed!", NULL);
  }
}

namespace SmartMet
{
namespace
{
// ----------------------------------------------------------------------
/*!
 * \brief Get points on an isocircle
 */
// ----------------------------------------------------------------------

std::pair<checkedVector<NFmiPoint>, std::vector<double>> get_isocircle_points(
    double lon1, double lat1, double lon2, double lat2, std::size_t steps)
{
  try
  {
    // Sanity checks
    if (lon1 == lon2 && lat1 == lat2)
      throw Spine::Exception(BCP, "Ill-defined isocircle: start and end points are equal");

    if (std::abs(lon1 - lon2) == 180 && std::abs(lat1 - (90 - lat2)) == 90)
      throw Spine::Exception(BCP, "Ill-defined isocircle: points at opposing ends of the earth");

    if (steps < 1 || steps > 10000)
      throw Spine::Exception(BCP,
                             "Number of points on isocircle must be 1-10000, not " +
                                 boost::lexical_cast<std::string>(steps));

    // Calculate bearing and distance to be travelled

    NFmiLocation startpoint(lon1, lat1);
    NFmiLocation endpoint(lon2, lat2);

    double bearing = startpoint.Direction(endpoint);
    double distance = startpoint.Distance(endpoint);

    checkedVector<NFmiPoint> coordinates;
    coordinates.push_back(startpoint.GetLocation());

    std::vector<double> distances;
    distances.push_back(0);

    for (std::size_t i = 1; i < steps; i++)
    {
      // Should this be fixed? Probably not - the coordinates should behave the same
      const bool pacific_view = false;
      double dist = i * distance / steps;
      auto loc = startpoint.GetLocation(bearing, dist, pacific_view);
      coordinates.push_back(loc.GetLocation());
      distances.push_back(dist / 1000.0);
    }
    coordinates.push_back(endpoint.GetLocation());
    distances.push_back(distance / 1000.0);

    return std::make_pair(coordinates, distances);
  }
  catch (...)
  {
    throw Spine::Exception(BCP, "Operation failed!", NULL);
  }
}

}  // namespace

namespace Engine
{
namespace Contour
{
// ----------------------------------------------------------------------
/*!
 * \brief Convert a GEOS geometry to OGR form
 */
// ----------------------------------------------------------------------

OGRGeometryPtr geos_to_ogr(const GeometryPtr &theGeom, OGRSpatialReference *theSR)
{
  try
  {
    if (!theGeom)
      return OGRGeometryPtr();

    // Convert to WKB

    std::ostringstream out;
    geos::io::WKBWriter writer;
    writer.write(*theGeom, out);
    const auto &wkb = out.str();

    // Read to OGR using as hideous casts as is required

    unsigned char *cwkb = reinterpret_cast<unsigned char *>(const_cast<char *>(wkb.c_str()));

    OGRGeometry *ogeom;
    OGRErr err = OGRGeometryFactory::createFromWkb(cwkb, theSR, &ogeom);

    if (err != OGRERR_NONE)
    {
      // we assume createFromWkb does not destroy the spatial reference on failure
      delete theSR;
      throw Spine::Exception(BCP, "Failed to convert contoured WKB to OGRGeometry");
    }

    if (theSR != NULL)
      theSR->Dereference();  // Spatial references are reference counted, we must relinquish the
                             // local
                             // copy here, otherwise this leaks

    // Return the result
    return OGRGeometryPtr(ogeom);
  }
  catch (...)
  {
    throw Spine::Exception(BCP, "Operation failed!", NULL);
  }
}

// ----------------------------------------------------------------------
/*!
 * \brief Implementation constructor
 */
// ----------------------------------------------------------------------

Engine::Impl::Impl(const std::string &theFileName)
    : itsConfigFile(theFileName), itsGeomFactory(new geos::geom::GeometryFactory())
{
}

// ----------------------------------------------------------------------
/*!
 * \brief Initialization is done in a separate thread
 */
// ----------------------------------------------------------------------

void Engine::Impl::init()
{
  try
  {
    itsConfig.reset(new Config(itsConfigFile));
    itsContourCache.resize(itsConfig->getMaxContourCacheSize());
  }
  catch (...)
  {
    throw Spine::Exception(BCP, "Operation failed!", NULL);
  }
}

// ----------------------------------------------------------------------
/*!
 * \brief Get cache usage report
 */
// ----------------------------------------------------------------------

CacheReportingStruct Engine::Impl::getCacheSizes()
{
  try
  {
    CacheReportingStruct ret;

    ret.contour_cache_max_size = itsContourCache.maxSize();
    ret.contour_cache_size = itsContourCache.size();

    return ret;
  }
  catch (...)
  {
    throw Spine::Exception(BCP, "Operation failed!", NULL);
  }
}

// ----------------------------------------------------------------------
/*!
 * \brief Empty the cache
 */
// ----------------------------------------------------------------------

void Engine::Impl::clearCache()
{
  itsContourCache.clear();
}

// ----------------------------------------------------------------------
/*!
 * \brief Low level contouring utility
 *
 * Note: We intentionally copy the values to be able to filter the data.
 *       We cache the data values, hence modifying the input is not OK.
 */
// ----------------------------------------------------------------------

GeometryPtr Engine::Impl::internal_isoline(const DataMatrixAdapter &data,
                                           const MyHints &hints,
                                           bool worldwrap,
                                           double isovalue,
                                           Interpolation interpolation)
{
  // Should support multiple builders with different SRIDs
  Tron::FmiBuilder builder(itsGeomFactory);

  // isoline
  switch (interpolation)
  {
    case Linear:
    {
      MyLinearContourer::line(builder, data, isovalue, worldwrap, hints);
      break;
    }
    case LogLinear:
    {
      MyLogLinearContourer::line(builder, data, isovalue, worldwrap, hints);
      break;
    }
    case Nearest:
    {
      throw std::runtime_error("Isolines not supported with nearest neighbour interpolation");
    }
    case Discrete:
    {
      throw std::runtime_error("Isolines not supported with discrete interpolation");
    }
  }
  return builder.result();
}

// ----------------------------------------------------------------------
/*!
 * \brief Low level contouring utility
 *
 * Note: We intentionally copy the values to be able to filter the data.
 *       We cache the data values, hence modifying the input is not OK.
 */
// ----------------------------------------------------------------------

GeometryPtr Engine::Impl::internal_isoband(const DataMatrixAdapter &data,
                                           const MyHints &hints,
                                           bool worldwrap,
                                           const boost::optional<double> &lolimit,
                                           const boost::optional<double> &hilimit,
                                           Interpolation interpolation)
{
  // Should support multiple builders with different SRIDs
  Tron::FmiBuilder builder(itsGeomFactory);

  double lo = kFloatMissing, hi = kFloatMissing;
  if (lolimit)
    lo = *lolimit;
  if (hilimit)
    hi = *hilimit;

  switch (interpolation)
  {
    case Linear:
    {
      MyLinearContourer::fill(builder, data, lo, hi, worldwrap, hints);
      break;
    }
    case LogLinear:
    {
      MyLogLinearContourer::fill(builder, data, lo, hi, worldwrap, hints);
      break;
    }
    case Nearest:
    {
      MyNearestContourer::fill(builder, data, lo, hi, worldwrap, hints);
      break;
    }
    case Discrete:
    {
      MyDiscreteContourer::fill(builder, data, lo, hi, worldwrap, hints);
      break;
    }
  }

  return builder.result();
}

// ----------------------------------------------------------------------
/*!
 * \brief Contour producing vector of OGR geometries
 */
// ----------------------------------------------------------------------

std::vector<OGRGeometryPtr> Engine::Impl::contour(std::size_t theQhash,
                                                  const std::string &theQAreaWKT,
                                                  const NFmiDataMatrix<float> &theMatrix,
                                                  const CoordinatesPtr theCoordinates,
                                                  const Options &theOptions,
                                                  bool worldwrap,
                                                  OGRSpatialReference *theSR)
{
  try
  {
    // Safety check, the result will be completely empty. Returning empty
    // OGRGeometryPtr values is not what for example the WFS plugin expects,
    // it expects empty geometries instead. It can however handle an empty
    // vector as the return value.

    std::size_t nx = theMatrix.NX();
    std::size_t ny = theMatrix.NY();

    if (nx == 0 || ny == 0)
      return {};

    // The hash for the result
    auto common_hash = theQhash;
    boost::hash_combine(common_hash, theOptions.filtered_data_hash_value());
    if (theSR != nullptr)
      boost::hash_combine(common_hash, *theSR);

    // results first include isolines, then isobands
    auto nresults = theOptions.isovalues.size() + theOptions.limits.size();

    // The results
    std::vector<OGRGeometryPtr> retval(nresults);

    // Make a copy of input data to enable filtering.
    // TODO: Use lazy initialization as in data and hints

    NFmiDataMatrix<float> values = theMatrix;

    // We generate helper data only if needed

    std::unique_ptr<DataMatrixAdapter> data;
    std::unique_ptr<MyHints> hints;

    // Prepare the SR from querydata only once

    std::unique_ptr<OGRSpatialReference> querydata_sr;

    if (theSR == nullptr)
    {
      querydata_sr.reset(new OGRSpatialReference);
      OGRErr err = querydata_sr->SetFromUserInput(theQAreaWKT.c_str());
      if (err != OGRERR_NONE)
        throw std::runtime_error("GDAL does not understand this FMI WKT: " + theQAreaWKT);
    }

    // Parallel processing of the contours

    struct Args
    {
      std::size_t i;
      std::size_t hash;
    };

    const auto contourer =
        [&retval, this, &theOptions, &theSR, &querydata_sr, &worldwrap, &data, &hints](Args &args) {
          try
          {
            // Calculate GEOS geometry result
            GeometryPtr geom;
            if (args.i < theOptions.isovalues.size())
            {
              geom = internal_isoline(
                  *data, *hints, worldwrap, theOptions.isovalues[args.i], theOptions.interpolation);
            }
            else
            {
              auto iband = args.i - theOptions.isovalues.size();
              geom = internal_isoband(*data,
                                      *hints,
                                      worldwrap,
                                      theOptions.limits[iband].lolimit,
                                      theOptions.limits[iband].hilimit,
                                      theOptions.interpolation);
            }

            // Cache as OGRGeometry

            if (!geom)
            {
              // We want to cache empty results too to save speed!
              auto ret = OGRGeometryPtr();
              itsContourCache.insert(args.hash, ret);
              retval[args.i] = ret;
              return;
            }

            // Generate spatial reference. This needs to be done with new,
            // createFromWkb will take ownership of the pointer passed to it.

            OGRSpatialReference *sr = nullptr;
            if (theSR != nullptr)
              sr = theSR->Clone();
            else
              sr = querydata_sr->Clone();

            // Convert to OGR object and cache the result

            auto ret = geos_to_ogr(geom, sr);
            retval[args.i] = ret;
            itsContourCache.insert(args.hash, ret);
          }
          catch (...)
          {
            // Cannot let an exception cause termination. We'll let the user get an empty result
            // instead.
            Spine::Exception ex(BCP, "Contouring failed");
            ex.printError();
          }
        };

    // Do not use all cores until better balancing & locking is implemented
    const auto max_concurrency = std::thread::hardware_concurrency();
    const auto concurrency = std::max(1u, max_concurrency / 4);

    Fmi::WorkQueue<Args> workqueue(contourer, concurrency);

    for (auto icontour = 0ul; icontour < nresults; icontour++)
    {
      // Establish cache hash
      auto hash = common_hash;
      if (icontour < theOptions.isovalues.size())
      {
        boost::hash_combine(hash, boost::hash_value("isoline"));
        boost::hash_combine(hash, boost::hash_value(theOptions.isovalues[icontour]));
      }
      else
      {
        auto i = icontour - theOptions.isovalues.size();
        boost::hash_combine(hash, boost::hash_value("isoband"));

        // TODO: Use generic code as in wms/Hash.h
        if (theOptions.limits[i].lolimit)
          boost::hash_combine(hash, boost::hash_value(*theOptions.limits[i].lolimit));
        else
          boost::hash_combine(hash, boost::hash_value(false));

        if (theOptions.limits[i].hilimit)
          boost::hash_combine(hash, boost::hash_value(*theOptions.limits[i].hilimit));
        else
          boost::hash_combine(hash, boost::hash_value(false));
      }

      auto opt_geom = itsContourCache.find(hash);
      if (opt_geom)
      {
        retval[icontour] = *opt_geom;
        continue;
      }

      // Perform unit conversion, smoothing and building hints only once.
      // This will be done only if some contour wasn't in the cache.

      if (!data)
      {
        // Perform unit conversion

        if (theOptions.hasTransformation())
        {
          double a = (theOptions.multiplier ? *theOptions.multiplier : 1.0);
          double b = (theOptions.offset ? *theOptions.offset : 0.0);

          std::size_t nx = values.NX();
          std::size_t ny = values.NY();

          for (std::size_t j = 0; j < ny; j++)
            for (std::size_t i = 0; i < nx; i++)
              if (values[i][j] != kFloatMissing)
                values[i][j] = a * values[i][j] + b;
        }

        // Adapter for contouring
        data.reset(new DataMatrixAdapter(values, *theCoordinates));

        if (theOptions.filter_size || theOptions.filter_degree)
        {
          size_t size = (theOptions.filter_size ? *theOptions.filter_size : 1);
          size_t degree = (theOptions.filter_degree ? *theOptions.filter_degree : 1);
          Tron::SavitzkyGolay2D::smooth(*data, size, degree);
        }

        // Search tree for extremum values
        hints.reset(new MyHints(*data));
      }

      // Calculate GEOS geometry result in a separate thread
      Args args{icontour, hash};
      workqueue(args);
    }

    workqueue.join_all();

    return retval;
  }
  catch (...)
  {
    throw Spine::Exception(BCP, "Operation failed!", NULL);
  }
}

// ----------------------------------------------------------------------
/*!
 * \brief Produce a cross section contur
 */
// ----------------------------------------------------------------------

std::vector<OGRGeometryPtr> Engine::Impl::crossection(
    boost::shared_ptr<NFmiFastQueryInfo> theQInfo,
    const Options &theOptions,
    const boost::optional<Spine::Parameter> &theZParameter,
    double theLon1,
    double theLat1,
    double theLon2,
    double theLat2,
    std::size_t theSteps)
{
  try
  {
    // Trivial option checks

    if (theOptions.level)
      throw Spine::Exception(BCP, "Using the level option is meaningless for cross-sections");
    if (theOptions.filter_size)
      throw Spine::Exception(BCP, "Using the filter_size option is meaningless for cross-sections");
    if (theOptions.filter_degree)
      throw Spine::Exception(BCP,
                             "Using the filter_degree option is meaningless for cross-sections");

    // Verify height parameter is available

    unsigned long zparam = 0;

    if (theZParameter)
    {
      if (!theQInfo->Param(theZParameter->number()))
      {
        std::cerr << "CSection: ZParameter " << theZParameter->name() << " is not available"
                  << std::endl;
        // TODO: Give good error message
        return {};
      }
      zparam = theQInfo->ParamIndex();
    }

    // Select the data

    if (!theQInfo->Param(theOptions.parameter.number()))
    {
      std::cerr << "CSection: Parameter " << theOptions.parameter.name() << " is not available"
                << std::endl;
      // TODO: Give good error message
      return {};
    }

    unsigned long param = theQInfo->ParamIndex();

    // Calculate the isocircle points

    const auto isocircle = get_isocircle_points(theLon1, theLat1, theLon2, theLat2, theSteps);
    const auto &coordinates = isocircle.first;
    const auto &distances = isocircle.second;

    // Get the cross section

    // TODO: Q API MUST BE CLEANED UP!

    if (!theQInfo->IsInside(theOptions.time))
    {
      // TODO: Give good error message
      return {};
    }

    if (theQInfo->SizeLevels() == 1)
    {
      // TODO: Give good error message
      return {};
    }

    std::size_t nx = coordinates.size();
    std::size_t ny = theQInfo->SizeLevels();

    NFmiDataMatrix<float> values(nx, ny, kFloatMissing);
    NFmiDataMatrix<NFmiPoint> coords(nx, ny, NFmiPoint());

    bool has_some_valid_levelvalues = false;

    std::size_t j = 0;
    for (theQInfo->ResetLevel(); theQInfo->NextLevel(); ++j)
    {
      float levelvalue = theQInfo->Level()->LevelValue();
      for (std::size_t i = 0; i < coordinates.size(); i++)
      {
        const int max_minutes = 360;
        values[i][j] = theQInfo->InterpolatedValue(coordinates[i], theOptions.time, max_minutes);

        if (theZParameter)
        {
          theQInfo->ParamIndex(zparam);
          levelvalue = theQInfo->InterpolatedValue(coordinates[i], theOptions.time, max_minutes);
          theQInfo->ParamIndex(param);
          // TODO: Handle missing values (how??)
        }

        if (levelvalue == kFloatMissing)
          levelvalue = std::numeric_limits<float>::quiet_NaN();

        if (!std::isnan(levelvalue))
          has_some_valid_levelvalues = true;

        coords[i][j] = NFmiPoint(distances[i], levelvalue);
      }
    }

    if (!has_some_valid_levelvalues)
      throw std::runtime_error("The heights for the cross-section are missing");

    // Make sure the y-coordinate is increasing, it may be decreasing for pressure
    // level or model level data.

    if (coords[0][1] < coords[0][0])
    {
      for (std::size_t i = 0; i < nx; i++)
      {
        std::size_t j1 = 0;
        std::size_t j2 = ny - 1;
        while (j1 < j2)
        {
          std::swap(values[i][j1], values[i][j2]);
          std::swap(coords[i][j1], coords[i][j2]);
          ++j1;
          --j2;
        }
      }
    }

    DataMatrixAdapter data(values, coords);
    MyHints hints(data);

    // results first include isolines, then isobands
    auto nresults = theOptions.isovalues.size() + theOptions.limits.size();

    // The results
    std::vector<OGRGeometryPtr> retval(nresults);

    for (auto icontour = 0ul; icontour < nresults; icontour++)
    {
      // Is it an isovalue request?
      bool isovaluerequest = (icontour < theOptions.isovalues.size());

      bool worldwrap = false;

      // Calculate GEOS geometry result
      GeometryPtr geom;
      if (isovaluerequest)
      {
        geom = internal_isoline(
            data, hints, worldwrap, theOptions.isovalues[icontour], theOptions.interpolation);
      }
      else
      {
        auto i = icontour - theOptions.isovalues.size();
        geom = internal_isoband(data,
                                hints,
                                worldwrap,
                                theOptions.limits[i].lolimit,
                                theOptions.limits[i].hilimit,
                                theOptions.interpolation);
      }

      OGRSpatialReference *sr = NULL;
      retval[icontour] = geos_to_ogr(geom, sr);
    }

    return retval;
  }
  catch (...)
  {
    throw Spine::Exception(BCP, "Operation failed!", NULL);
  }
}

}  // namespace Contour
}  // namespace Engine
}  // namespace SmartMet
