// ======================================================================
/*!
 * \brief Implementation details for the Contour engine
 */
// ======================================================================

#include "Engine.h"
#include "Config.h"
#include "DataMatrixAdapter.h"
#include "Engine.h"
#include "GeosTools.h"
#include "Options.h"
#include <boost/functional/hash.hpp>
#include <geos/geom/Geometry.h>
#include <geos/geom/GeometryFactory.h>
#include <geos/io/WKBWriter.h>
#include <geos/version.h>
#include <gis/CoordinateMatrixAnalysis.h>
#include <gis/OGR.h>
#include <macgyver/Cache.h>
#include <macgyver/StringConversion.h>
#include <macgyver/WorkQueue.h>
#include <spine/Exception.h>
#include <tron/FmiBuilder.h>
#include <tron/SavitzkyGolay2D.h>
#include <tron/Tron.h>
#include <cmath>
#include <cpl_conv.h>  // For configuring GDAL
#include <limits>
#include <memory>
#include <ogr_core.h>
#include <ogr_geometry.h>
#include <ogr_spatialref.h>

using GeometryPtr = std::shared_ptr<geos::geom::Geometry>;

// Contourers

using MyTraits = Tron::Traits<double, double, Tron::InfMissing>;
using MyBuilder = Tron::FmiBuilder;
using MyHints = Tron::Hints<DataMatrixAdapter, MyTraits>;

using MyLinearContourer =
    Tron::Contourer<DataMatrixAdapter, MyBuilder, MyTraits, Tron::LinearInterpolation>;

using MyLogLinearContourer =
    Tron::Contourer<DataMatrixAdapter, MyBuilder, MyTraits, Tron::LogLinearInterpolation>;

using MyNearestContourer =
    Tron::Contourer<DataMatrixAdapter, MyBuilder, MyTraits, Tron::NearestNeighbourInterpolation>;

using MyDiscreteContourer =
    Tron::Contourer<DataMatrixAdapter, MyBuilder, MyTraits, Tron::DiscreteInterpolation>;

namespace SmartMet
{
namespace Engine
{
namespace Contour
{
class Engine::Impl
{
 public:
  Impl(const std::string &theFileName);
  Impl() = delete;

  void init();

  // Produce vector of OGR contours for the given spatial reference

  std::vector<OGRGeometryPtr> contour(std::size_t theQhash,
                                      const Fmi::SpatialReference &theDataCRS,
                                      const Fmi::SpatialReference &theOutputCRS,
                                      const NFmiDataMatrix<float> &theMatrix,
                                      const CoordinatesPtr theCoordinates,
                                      const Options &theOptions,
                                      bool worldwrap) const;

  // Produce an OGR crossection for the given data
  std::vector<OGRGeometryPtr> crossection(NFmiFastQueryInfo &theQInfo,
                                          const Options &theOptions,
                                          const boost::optional<Spine::Parameter> &theZParameter,
                                          double theLon1,
                                          double theLat1,
                                          double theLon2,
                                          double theLat2,
                                          std::size_t theSteps);

  CacheReportingStruct getCacheSizes();
  void clearCache();

 private:
  std::string itsConfigFile;
  std::unique_ptr<Config> itsConfig;  // ptr for delayed initialization

  geos::geom::GeometryFactory::Ptr itsGeomFactory;

  // Cached contours

  using GeometryCache = Fmi::Cache::Cache<std::size_t, OGRGeometryPtr>;
  mutable GeometryCache itsContourCache;

  // Cached information on the handedness of grids
  using HandednessCache = Fmi::Cache::Cache<std::size_t, bool>;
  mutable HandednessCache itsHandednessCache;

  bool needs_flipping(const Fmi::CoordinateMatrix &coords,
                      std::size_t datahash,
                      OGRSpatialReference *srs) const;

  GeometryPtr internal_isoline(const DataMatrixAdapter &data,
                               const MyHints &hints,
                               bool worldwrap,
                               double isovalue,
                               Interpolation interpolation);

  GeometryPtr internal_isoband(const DataMatrixAdapter &data,
                               const MyHints &hints,
                               bool worldwrap,
                               const boost::optional<double> &lolimit,
                               const boost::optional<double> &hilimit,
                               Interpolation interpolation);
};

}  // namespace Contour
}  // namespace Engine
}  // namespace SmartMet

// IMPLEMENTATION DETAILS

namespace
{
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
    CPLFree(wkt);
    boost::hash<std::string> hasher;
    return hasher(tmp);
  }
  catch (...)
  {
    throw SmartMet::Spine::Exception::Trace(BCP, "Operation failed!");
  }
}

// ----------------------------------------------------------------------
/*!
 * \brief Utility function for extrapolation
 */
// ----------------------------------------------------------------------

class Extrapolation
{
 public:
  bool ok() { return count > 0; }
  float result() const { return sum / count; }
  void operator()(float value)
  {
    if (!std::isnan(value))
    {
      sum += value;
      ++count;
    }
  }

 private:
  float sum = 0;
  int count = 0;
};

// ----------------------------------------------------------------------
/*!
 * \brief Extrapolate data.
 *
 * This is typically used to expand the data a bit when it covers for example
 * only land areas, but visualization would leave gaps between the data and the
 * shoreline.
 /*/

void extrapolate(NFmiDataMatrix<float> &theValues, int theAmount)
{
  if (theAmount <= 0)
    return;

  auto tmp = theValues;
  for (int i = 0; i < theAmount; i++)
  {
    for (unsigned int j = 0; j < theValues.NY(); j++)
      for (unsigned int i = 0; i < theValues.NX(); i++)
      {
        if (std::isnan(tmp[i][j]))
        {
          Extrapolation ext;
          ext(theValues.At(i - 1, j, std::numeric_limits<float>::quiet_NaN()));
          ext(theValues.At(i + 1, j, std::numeric_limits<float>::quiet_NaN()));
          ext(theValues.At(i, j - 1, std::numeric_limits<float>::quiet_NaN()));
          ext(theValues.At(i, j + 1, std::numeric_limits<float>::quiet_NaN()));
          if (!ext.ok())
          {
            ext(theValues.At(i - 1, j - 1, std::numeric_limits<float>::quiet_NaN()));
            ext(theValues.At(i - 1, j + 1, std::numeric_limits<float>::quiet_NaN()));
            ext(theValues.At(i + 1, j - 1, std::numeric_limits<float>::quiet_NaN()));
            ext(theValues.At(i + 1, j + 1, std::numeric_limits<float>::quiet_NaN()));
          }
          if (ext.ok())
            tmp[i][j] = ext.result();
        }
      }

    std::swap(tmp, theValues);
  }
}

}  // anonymous namespace

namespace SmartMet
{
namespace
{
// ----------------------------------------------------------------------
/*!
 * \brief Get points on an isocircle
 */
// ----------------------------------------------------------------------

std::pair<std::vector<NFmiPoint>, std::vector<double>> get_isocircle_points(
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

    std::vector<NFmiPoint> coordinates;
    coordinates.push_back(startpoint.GetLocation());

    std::vector<double> distances;
    distances.push_back(0);

    for (std::size_t i = 1; i < steps; i++)
    {
      // Should this be fixed? Probably not - the coordinates should behave the same
      double dist = i * distance / steps;
      auto loc = startpoint.GetLocation(bearing, dist);
      coordinates.push_back(loc.GetLocation());
      distances.push_back(dist / 1000.0);
    }
    coordinates.push_back(endpoint.GetLocation());
    distances.push_back(distance / 1000.0);

    return std::make_pair(coordinates, distances);
  }
  catch (...)
  {
    throw Spine::Exception::Trace(BCP, "Operation failed!");
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

    if (theSR != nullptr)
      theSR->Dereference();  // Spatial references are reference counted, we must relinquish the
                             // local
                             // copy here, otherwise this leaks

    // Return the result
    return OGRGeometryPtr(ogeom);
  }
  catch (...)
  {
    throw Spine::Exception::Trace(BCP, "Operation failed!");
  }
}

// ----------------------------------------------------------------------
/*!
 * \brief Implementation constructor
 */
// ----------------------------------------------------------------------

Engine::Impl::Impl(const std::string &theFileName)
    : itsConfigFile(theFileName), itsGeomFactory(geos::geom::GeometryFactory::create())
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

    // 10,000 booleans doesn't take much space
    itsHandednessCache.resize(10000);
  }
  catch (...)
  {
    throw Spine::Exception::Trace(BCP, "Operation failed!");
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
    throw Spine::Exception::Trace(BCP, "Operation failed!");
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
  Tron::FmiBuilder builder(*itsGeomFactory);

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
  Tron::FmiBuilder builder(*itsGeomFactory);

  double lo = -std::numeric_limits<double>::infinity();
  double hi = +std::numeric_limits<double>::infinity();

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
 * \brief Returns true if the handedness of the data is incorrect
 */
// ----------------------------------------------------------------------

bool Engine::Impl::needs_flipping(const Fmi::CoordinateMatrix &coords,
                                  std::size_t datahash,
                                  OGRSpatialReference *crs) const
{
  // Combined hash value includes data details and output projection
  std::size_t hash = datahash;
  if (crs != nullptr)
    boost::hash_combine(hash, boost::hash_value(Fmi::OGR::exportToWkt(*crs)));

  // See if the result has already been cached

  auto cached = itsHandednessCache.find(hash);
  if (cached)
    return *cached;

  // Calculate handedness

  bool handedness = false;

  // Cache the result

  itsHandednessCache.insert(hash, handedness);
  return handedness;
}

// ----------------------------------------------------------------------
/*!
 * \brief Change all kFloatMissing values to NaN
 */
// ----------------------------------------------------------------------

void set_missing_to_nan(NFmiDataMatrix<float> &theValues)
{
  const std::size_t nx = theValues.NX();
  const std::size_t ny = theValues.NY();
  if (nx == 0 || ny == 0)
    return;

  const auto nan = std::numeric_limits<float>::quiet_NaN();

  // Unfortunately NFmiDataMatrix is a vector of vectors, memory
  // access patterns are not optimal

  for (std::size_t i = 0; i < nx; i++)
  {
    auto &tmp = theValues[i];
    for (std::size_t j = 0; j < ny; j++)
    {
      if (tmp[j] == kFloatMissing)
        tmp[j] = nan;
    }
  }
}

// ----------------------------------------------------------------------
/*!
 * \brief Calculate cache keys for the contours to be calculated
 */
// ----------------------------------------------------------------------

std::vector<std::size_t> contour_cache_keys(std::size_t theQhash,
                                            const Fmi::SpatialReference &theOutputCRS,
                                            const Options &theOptions)
{
  // The hash for the result
  auto common_hash = theQhash;
  boost::hash_combine(common_hash, theOptions.filtered_data_hash_value());
  boost::hash_combine(common_hash, theOutputCRS.hashValue());

  // results first include isolines, then isobands
  const auto nresults = theOptions.isovalues.size() + theOptions.limits.size();

  std::vector<std::size_t> hash_values(nresults, 0);

  const auto isoline_hash = boost::hash_value("isoline");
  const auto isoband_hash = boost::hash_value("isoband");

  for (auto icontour = 0ul; icontour < nresults; icontour++)
  {
    auto hash = common_hash;
    if (icontour < theOptions.isovalues.size())
    {
      boost::hash_combine(hash, isoline_hash);
      boost::hash_combine(hash, boost::hash_value(theOptions.isovalues[icontour]));
    }
    else
    {
      boost::hash_combine(hash, isoband_hash);

      auto i = icontour - theOptions.isovalues.size();

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

    hash_values[icontour] = hash;
  }
  return hash_values;
}

// ----------------------------------------------------------------------
/*!
 * \brief Contour producing vector of OGR geometries
 */
// ----------------------------------------------------------------------

std::vector<OGRGeometryPtr> Engine::Impl::contour(std::size_t theQhash,
                                                  const Fmi::SpatialReference &theDataCRS,
                                                  const Fmi::SpatialReference &theOutputCRS,
                                                  const NFmiDataMatrix<float> &theMatrix,
                                                  const CoordinatesPtr theCoordinates,
                                                  const Options &theOptions,
                                                  bool worldwrap) const
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

    // results first include isolines, then isobands
    auto nresults = theOptions.isovalues.size() + theOptions.limits.size();

    std::vector<OGRGeometryPtr> retval(nresults);

    // Calculate the cache keys for all the contours. This enables us to test
    // if everything is cached already without doing any data processing.

    std::vector<std::size_t> contour_hash_values =
        contour_cache_keys(theQhash, theOutputCRS, theOptions);

    // Make a copy of input data to enable filtering.
    // TODO: Use lazy initialization as in data and hints

    NFmiDataMatrix<float> values = theMatrix;
    CoordinatesPtr coords = theCoordinates;

    // Make sure the values are NaN instead of kFloatMissing
    set_missing_to_nan(values);

    if (needs_flipping(*coords, theQhash, theDataCRS))
    {
      for (std::size_t j1 = 0; j1 < ny / 2; j1++)
      {
        std::size_t j2 = ny - j1 - 1;

        for (std::size_t i = 0; i < nx; i++)
          std::swap(values[i][j1], values[i][j2]);
      }

      // we copy the coordinates only when we have to flip them
      coords.reset(new Fmi::CoordinateMatrix(*coords));

      for (std::size_t j1 = 0; j1 < ny / 2; j1++)
      {
        std::size_t j2 = ny - j1 - 1;

        for (std::size_t i = 0; i < nx; i++)
        {
          NFmiPoint pt1 = (*coords)(i, j1);
          NFmiPoint pt2 = (*coords)(i, j2);
          coords->set(i, j1, pt2);
          coords->set(i, j2, pt1);
        }
      }
    }

    // We generate helper data only if needed

    std::unique_ptr<DataMatrixAdapter> data;
    std::unique_ptr<MyHints> hints;

    // Parallel processing of the contours

    struct Args
    {
      std::size_t i;
      std::size_t hash;
    };

    const auto contourer =
        [&retval, this, &theOptions, &theDataCRS, &theOutputCRS, &worldwrap, &data, &hints](
            Args &args) {
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

            // Convert to OGR object and cache the result

            auto ret = geos_to_ogr(geom, theOutputCRS.get()->Clone());

            // Despeckle even closed isolines (pressure curves)
            if (theOptions.minarea)
              ret.reset(Fmi::OGR::despeckle(*ret, *theOptions.minarea));

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
    const auto concurrency = std::max(1u, max_concurrency / 8);

    Fmi::WorkQueue<Args> workqueue(contourer, concurrency);

    for (auto icontour = 0ul; icontour < nresults; icontour++)
    {
      // Establish cache hash
      auto hash = contour_hash_values[icontour];

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
              values[i][j] = a * values[i][j] + b;
        }

        // Extrapolate if requested
        extrapolate(values, theOptions.extrapolation);

        // Adapter for contouring
        data.reset(
            new DataMatrixAdapter(values, *coords, Fmi::analysis(*coords).m_invalid));  // SHIT

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
    throw Spine::Exception::Trace(BCP, "Operation failed!");
  }
}

// ----------------------------------------------------------------------
/*!
 * \brief Produce a cross section contur
 */
// ----------------------------------------------------------------------

std::vector<OGRGeometryPtr> Engine::Impl::crossection(
    NFmiFastQueryInfo &theQInfo,
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
      if (!theQInfo.Param(theZParameter->number()))
      {
        std::cerr << "CSection: ZParameter " << theZParameter->name() << " is not available"
                  << std::endl;
        // TODO: Give good error message
        return {};
      }
      zparam = theQInfo.ParamIndex();
    }

    // Select the data

    if (!theQInfo.Param(theOptions.parameter.number()))
    {
      std::cerr << "CSection: Parameter " << theOptions.parameter.name() << " is not available"
                << std::endl;
      // TODO: Give good error message
      return {};
    }

    unsigned long param = theQInfo.ParamIndex();

    // Calculate the isocircle points

    const auto isocircle = get_isocircle_points(theLon1, theLat1, theLon2, theLat2, theSteps);
    const auto &coordinates = isocircle.first;
    const auto &distances = isocircle.second;

    // Get the cross section

    // TODO: Q API MUST BE CLEANED UP!

    if (!theQInfo.IsInside(theOptions.time))
    {
      // TODO: Give good error message
      return {};
    }

    if (theQInfo.SizeLevels() == 1)
    {
      // TODO: Give good error message
      return {};
    }

    std::size_t nx = coordinates.size();
    std::size_t ny = theQInfo.SizeLevels();

    NFmiDataMatrix<float> values(nx, ny, std::numeric_limits<float>::quiet_NaN());
    Fmi::CoordinateMatrix coords(nx, ny);

    bool has_some_valid_levelvalues = false;

    std::size_t j = 0;
    for (theQInfo.ResetLevel(); theQInfo.NextLevel(); ++j)
    {
      float levelvalue = theQInfo.Level()->LevelValue();
      for (std::size_t i = 0; i < coordinates.size(); i++)
      {
        const int max_minutes = 360;
        values[i][j] = theQInfo.InterpolatedValue(coordinates[i], theOptions.time, max_minutes);

        if (theZParameter)
        {
          theQInfo.ParamIndex(zparam);
          levelvalue = theQInfo.InterpolatedValue(coordinates[i], theOptions.time, max_minutes);
          theQInfo.ParamIndex(param);
          // TODO: Handle missing values (how??)
        }

        if (levelvalue == kFloatMissing)
          levelvalue = std::numeric_limits<float>::quiet_NaN();

        if (!std::isnan(levelvalue))
          has_some_valid_levelvalues = true;

        coords.set(i, j, distances[i], levelvalue);
      }
    }

    if (!has_some_valid_levelvalues)
      throw std::runtime_error("The heights for the cross-section are missing");

    // Make sure the y-coordinate is increasing, it may be decreasing for pressure
    // level or model level data.

    if (coords(0, 1) < coords(0, 0))
    {
      for (std::size_t i = 0; i < nx; i++)
      {
        std::size_t j1 = 0;
        std::size_t j2 = ny - 1;
        while (j1 < j2)
        {
          std::swap(values[i][j1], values[i][j2]);
          auto pt1 = coords(i, j1);
          auto pt2 = coords(i, j2);
          coords.set(i, j1, pt2.first, pt2.second);
          coords.set(i, j2, pt1.first, pt2.second);
          ++j1;
          --j2;
        }
      }
    }

    // Make sure the values are NaN instead of kFloatMissing
    set_missing_to_nan(values);

    // Don't bother analyzing the grid for cross sections, the Z coordinates
    // should always be fine.

    Fmi::BoolMatrix grid_is_fine(coords.width(), coords.height(), true);
    DataMatrixAdapter data(values, coords, grid_is_fine);
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

      OGRSpatialReference *sr = nullptr;
      retval[icontour] = geos_to_ogr(geom, sr);
    }

    return retval;
  }
  catch (...)
  {
    throw Spine::Exception::Trace(BCP, "Operation failed!");
  }
}

}  // namespace Contour
}  // namespace Engine
}  // namespace SmartMet

// ======================================================================
/*!
 * \brief Implementation of Contour::Engine
 */
// ======================================================================

namespace SmartMet
{
namespace Engine
{
namespace Contour
{
// ----------------------------------------------------------------------
/*!
 * \brief The only permitted constructor requires a configfile
 */
// ----------------------------------------------------------------------

Engine::Engine(const std::string &theFileName) : itsImpl(new Impl(theFileName)) {}
// ----------------------------------------------------------------------
/*!
 * \brief Initialize the engine
 */
// ----------------------------------------------------------------------
void Engine::init()
{
  try
  {
    itsImpl->init();

    // Discard projected points which fall outside the valid area for the spatial reference
    // For example eureffin is not valid in the full area of EC Europe data
    CPLSetConfigOption("OGR_ENABLE_PARTIAL_REPROJECTION", "YES");
  }
  catch (...)
  {
    throw Spine::Exception::Trace(BCP, "Operation failed!");
  }
}

// ----------------------------------------------------------------------
/*!
 * \brief Get cache usage report
 */
// ----------------------------------------------------------------------

CacheReportingStruct Engine::getCacheSizes() const
{
  try
  {
    return itsImpl->getCacheSizes();
  }
  catch (...)
  {
    throw Spine::Exception::Trace(BCP, "Operation failed!");
  }
}

// ----------------------------------------------------------------------
/*!
 * \brief Empty the caches. Mostly used while testing.
 */
// ----------------------------------------------------------------------

void Engine::clearCache() const
{
  itsImpl->clearCache();
}

// ----------------------------------------------------------------------
/*!
 * \brief Shutdown the engine
 */
// ----------------------------------------------------------------------

void Engine::shutdown()
{
  std::cout << "  -- Shutdown requested (contour)\n";
}
// ----------------------------------------------------------------------
/*!
 * \brief Contour
 */
// ----------------------------------------------------------------------

std::vector<OGRGeometryPtr> Engine::contour(std::size_t theQhash,
                                            const Fmi::SpatialReference &theDataCRS,
                                            const Fmi::SpatialReference &theOutputCRS,
                                            const NFmiDataMatrix<float> &theMatrix,
                                            const CoordinatesPtr theCoordinates,
                                            const Options &theOptions,
                                            bool worldwrap) const
{
  try
  {
    return itsImpl->contour(
        theQhash, theDataCRS, theOutputCRS, theMatrix, theCoordinates, theOptions, worldwrap);
  }
  catch (...)
  {
    throw Spine::Exception::Trace(BCP, "Operation failed!");
  }
}

// ----------------------------------------------------------------------
/*!
 * \brief Cross section
 */
// ----------------------------------------------------------------------

std::vector<OGRGeometryPtr> Engine::crossection(NFmiFastQueryInfo &theQInfo,
                                                const Options &theOptions,
                                                double theLon1,
                                                double theLat1,
                                                double theLon2,
                                                double theLat2,
                                                std::size_t theSteps) const
{
  try
  {
    boost::optional<Spine::Parameter> zparam;
    return itsImpl->crossection(
        theQInfo, theOptions, zparam, theLon1, theLat1, theLon2, theLat2, theSteps);
  }
  catch (...)
  {
    throw Spine::Exception::Trace(BCP, "Operation failed!");
  }
}

// ----------------------------------------------------------------------
/*!
 * \brief Cross section
 */
// ----------------------------------------------------------------------

std::vector<OGRGeometryPtr> Engine::crossection(NFmiFastQueryInfo &theQInfo,
                                                const Spine::Parameter &theZParameter,
                                                const Options &theOptions,
                                                double theLon1,
                                                double theLat1,
                                                double theLon2,
                                                double theLat2,
                                                std::size_t theSteps) const
{
  try
  {
    boost::optional<Spine::Parameter> zparam(theZParameter);
    return itsImpl->crossection(
        theQInfo, theOptions, zparam, theLon1, theLat1, theLon2, theLat2, theSteps);
  }
  catch (...)
  {
    throw Spine::Exception::Trace(BCP, "Operation failed!");
  }
}

}  // namespace Contour
}  // namespace Engine
}  // namespace SmartMet

// DYNAMIC MODULE CREATION TOOLS

extern "C" void *engine_class_creator(const char *configfile, void * /* user_data */)
{
  return new SmartMet::Engine::Contour::Engine(configfile);
}

extern "C" const char *engine_name()
{
  return "Contour";
}
// ======================================================================
