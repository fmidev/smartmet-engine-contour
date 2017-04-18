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
 * \brief Calculate a hash for the request
 */
// ----------------------------------------------------------------------

std::size_t contour_hash(std::size_t theQhash,
                         const Options &theOptions,
                         OGRSpatialReference *theSR)
{
  try
  {
    std::size_t seed = theQhash;
    // boost::hash_combine(seed, theQ);  // do not hash the shared pointer itself!
    boost::hash_combine(seed, theOptions);
    if (theSR != nullptr)
      boost::hash_combine(seed, *theSR);
    return seed;
  }
  catch (...)
  {
    throw Spine::Exception(BCP, "Operation failed!", NULL);
  }
}

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

GeometryPtr Engine::Impl::internal_contour(std::size_t datahash,
                                           const Options &options,
                                           const DataMatrixAdapter &data,
                                           const MyHints &hints,
                                           bool worldwrap)
{
  // Should support multiple builders with different SRIDs
  Tron::FmiBuilder builder(itsGeomFactory);

  if (!options.isovalues.empty())
  {
    auto isovalue = options.isovalues[0];

    // isoline
    switch (options.interpolation)
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
  }
  else
  {
    // isoband

    double lo = kFloatMissing, hi = kFloatMissing;
    Range isoband = options.limits[0];
    if (!!isoband.lolimit)
      lo = *isoband.lolimit;
    if (!!isoband.hilimit)
      hi = *isoband.hilimit;

    switch (options.interpolation)
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
  }

  return builder.result();
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
 * \brief Contour producing a GEOS geometry
 */
// ----------------------------------------------------------------------

GeometryPtr Engine::Impl::geosContour(std::size_t theQhash,
                                      const Options &theOptions,
                                      const DataMatrixAdapter &theData,
                                      const MyHints &theHints,
                                      bool worldwrap)
{
  try
  {
    // The hash for the final data depends on filtering etc options too

    auto datahash = theQhash;
    boost::hash_combine(datahash, theOptions.filtered_data_hash_value());

    return internal_contour(datahash, theOptions, theData, theHints, worldwrap);
  }
  catch (...)
  {
    throw Spine::Exception(BCP, "Operation failed!", NULL);
  }
}

// ----------------------------------------------------------------------
/*!
 * \brief Contour producing vector of OGR geometries
 */
// ----------------------------------------------------------------------

std::vector<OGRGeometryPtr> Engine::Impl::contour(std::size_t theQhash,
                                                  const std::string theQAreaWKT,
                                                  const NFmiDataMatrix<float> &theMatrix,
                                                  const CoordinatesPtr theCoordinates,
                                                  const Options &theOptions,
                                                  bool worldwrap,
                                                  OGRSpatialReference *theSR)
{
  std::vector<OGRGeometryPtr> retval;

  try
  {
    bool processIsoValues = !theOptions.isovalues.empty();
    unsigned int vectorSize =
        (processIsoValues ? theOptions.isovalues.size() : theOptions.limits.size());
    Options options = theOptions;

    std::unique_ptr<DataMatrixAdapter> data;
    std::unique_ptr<MyHints> hints;
    NFmiDataMatrix<float> values = theMatrix;

    for (unsigned int i = 0; i < vectorSize; i++)
    {
      if (processIsoValues)
      {
        options.isovalues.clear();
        options.isovalues.push_back(theOptions.isovalues[i]);
      }
      else
      {
        options.limits.clear();
        options.limits.push_back(theOptions.limits[i]);
      }

      // Try the cache first

      auto geom_hash = contour_hash(theQhash, options, theSR);

      auto opt_geom = itsContourCache.find(geom_hash);
      if (opt_geom)
      {
        retval.push_back(*opt_geom);
        continue;
      }

      std::size_t nx = theMatrix.NX();
      std::size_t ny = theMatrix.NY();

      if (nx == 0 || ny == 0)
        return retval;

      // Perform unit conversion, smoothing and building hints only once
      // and only if some response wasn't in the cache.

      if (!data)
      {
        // Perform unit conversion

        if (options.hasTransformation())
        {
          double a = (options.multiplier ? *options.multiplier : 1.0);
          double b = (options.offset ? *options.offset : 0.0);

          std::size_t nx = values.NX();
          std::size_t ny = values.NY();

          for (std::size_t j = 0; j < ny; j++)
            for (std::size_t i = 0; i < nx; i++)
              if (values[i][j] != kFloatMissing)
                values[i][j] = a * values[i][j] + b;
        }

        // Adapter for contouring
        data.reset(new DataMatrixAdapter(values, *theCoordinates));

        if (options.filter_size || options.filter_degree)
        {
          size_t size = (options.filter_size ? *options.filter_size : 1);
          size_t degree = (options.filter_degree ? *options.filter_degree : 1);
          Tron::SavitzkyGolay2D::smooth(*data, size, degree);
        }
        hints.reset(new MyHints(*data));
      }

      // Do the contouring
      auto geom = geosContour(theQhash, options, *data, *hints, worldwrap);

      if (!geom)
      {
        // We want to cache empty results too to save speed!
        auto ret = OGRGeometryPtr();
        itsContourCache.insert(geom_hash, ret);
        retval.push_back(ret);
        continue;
      }

      // Generate spatial reference. This needs to be done with new,
      // createFromWkb will take ownership of the pointer passed to it.

      OGRSpatialReference *sr;
      if (theSR != nullptr)
        sr = theSR->Clone();
      else
      {
        sr = new OGRSpatialReference;
        OGRErr err = sr->SetFromUserInput(theQAreaWKT.c_str());
        if (err != OGRERR_NONE)
        {
          delete sr;
          throw std::runtime_error("GDAL does not understand this FMI WKT: " + theQAreaWKT);
        }
      }

      // Convert to OGR object and cache the result

      auto ret = geos_to_ogr(geom, sr);
      retval.push_back(ret);
      itsContourCache.insert(geom_hash, ret);
    }
  }
  catch (...)
  {
    throw Spine::Exception(BCP, "Operation failed!", NULL);
  }

  return retval;
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
    std::vector<OGRGeometryPtr> ret;
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
        return ret;
      }
      zparam = theQInfo->ParamIndex();
    }

    // Select the data

    if (!theQInfo->Param(theOptions.parameter.number()))
    {
      std::cerr << "CSection: Parameter " << theOptions.parameter.name() << " is not available"
                << std::endl;
      // TODO: Give good error message
      return ret;
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
      return ret;
    }

    if (theQInfo->SizeLevels() == 1)
    {
      // TODO: Give good error message
      return ret;
    }

    std::size_t nx = coordinates.size();
    std::size_t ny = theQInfo->SizeLevels();

    NFmiDataMatrix<float> values(nx, ny, kFloatMissing);
    NFmiDataMatrix<NFmiPoint> coords(nx, ny, NFmiPoint());

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
        coords[i][j] = NFmiPoint(distances[i], levelvalue);
      }
    }

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

    bool processIsoValues = !theOptions.isovalues.empty();
    unsigned int vectorSize =
        (processIsoValues ? theOptions.isovalues.size() : theOptions.limits.size());
    Options options = theOptions;

    for (unsigned int i = 0; i < vectorSize; i++)
    {
      if (processIsoValues)
      {
        options.isovalues.clear();
        options.isovalues.push_back(theOptions.isovalues[i]);
      }
      else
      {
        options.limits.clear();
        options.limits.push_back(theOptions.limits[i]);
      }

      const bool worldwrap = false;  // until we encounter testable global data
      GeometryPtr geom = internal_contour(0, options, data, hints, worldwrap);
      OGRSpatialReference *sr = NULL;
      ret.push_back(geos_to_ogr(geom, sr));
    }

    return ret;
  }
  catch (...)
  {
    throw Spine::Exception(BCP, "Operation failed!", NULL);
  }
}

}  // namespace Contour
}  // namespace Engine
}  // namespace SmartMet
