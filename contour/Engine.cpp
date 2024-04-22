// ======================================================================
/*!
 * \brief Implementation details for the Contour engine
 */
// ======================================================================

#include "Engine.h"
#include "BaseGrid.h"
#include "Config.h"
#include "GeosTools.h"
#include "Options.h"
#include "PaddedGrid.h"
#include "ShiftedGrid.h"
#include <boost/optional.hpp>
#include <geos/geom/Geometry.h>
#include <geos/geom/GeometryFactory.h>
#include <geos/io/WKBWriter.h>
#include <geos/version.h>
#include <gis/CoordinateMatrixAnalysis.h>
#include <gis/EPSGInfo.h>
#include <gis/OGR.h>
#include <macgyver/Cache.h>
#include <macgyver/Exception.h>
#include <macgyver/Hash.h>
#include <macgyver/StringConversion.h>
#include <trax/Contour.h>
#include <trax/Geos.h>
#include <trax/IsobandLimits.h>
#include <trax/IsolineValues.h>
#include <trax/OGR.h>
#include <cmath>
#include <cpl_conv.h>  // For configuring GDAL
#include <future>
#include <limits>
#include <memory>
#include <ogr_core.h>
#include <ogr_geometry.h>
#include <ogr_spatialref.h>

using GeometryPtr = std::shared_ptr<geos::geom::Geometry>;

namespace SmartMet
{
namespace Engine
{
namespace Contour
{
// Contourers

class Engine::Impl
{
 public:
  explicit Impl(std::string theFileName);
  Impl() = delete;

  void init();

  // Produce vector of OGR contours for the given spatial reference

  std::vector<OGRGeometryPtr> contour(std::size_t theDataHash,
                                      const Fmi::SpatialReference &theOutputCRS,
                                      const NFmiDataMatrix<float> &theMatrix,
                                      const Fmi::CoordinateMatrix &theCoordinates,
                                      const Options &theOptions) const;

  std::vector<OGRGeometryPtr> contour(std::size_t theDataHash,
                                      const Fmi::SpatialReference &theOutputCRS,
                                      const NFmiDataMatrix<float> &theMatrix,
                                      const Fmi::CoordinateMatrix &theCoordinates,
                                      const Fmi::Box &theClipBox,
                                      const Options &theOptions) const;

  std::vector<OGRGeometryPtr> contour(std::size_t theDataHash,
                                      const Fmi::SpatialReference &theOutputCRS,
                                      const NFmiDataMatrix<float> &theMatrix,
                                      const Fmi::CoordinateMatrix &theCoordinates,
                                      const Fmi::Box &theClipBox,
                                      bool all_valid,
                                      Options theOptions) const;

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
  Fmi::Cache::CacheStatistics getCacheStats() const;

 private:
  std::string itsConfigFile;
  std::unique_ptr<Config> itsConfig;  // ptr for delayed initialization

  geos::geom::GeometryFactory::Ptr itsGeomFactory;

  // Cached contours

  using GeometryCache = Fmi::Cache::Cache<std::size_t, OGRGeometryPtr>;
  mutable GeometryCache itsContourCache;

  // Cached grid analyses

  using CoordinateAnalysisPtr = std::shared_ptr<Fmi::CoordinateAnalysis>;
  using AnalysisCache = Fmi::Cache::Cache<std::size_t, std::shared_future<CoordinateAnalysisPtr>>;
  mutable AnalysisCache itsAnalysisCache;

  CoordinateAnalysisPtr get_analysis(std::size_t theDataHash,
                                     const Fmi::CoordinateMatrix &theCoordinates,
                                     const Fmi::SpatialReference &theOutputCRS) const;

  GeometryPtr internal_isoline(const Trax::Grid &data,
                               double isovalue,
                               Trax::InterpolationType interpolation) const;

  GeometryPtr internal_isoband(const Trax::Grid &data,
                               const boost::optional<double> &lolimit,
                               const boost::optional<double> &hilimit,
                               Trax::InterpolationType interpolation) const;
};

}  // namespace Contour
}  // namespace Engine
}  // namespace SmartMet

// IMPLEMENTATION DETAILS

namespace
{
// ----------------------------------------------------------------------
/*!
 * \brief Utility class for extrapolation
 */
// ----------------------------------------------------------------------

class Extrapolation
{
 public:
  bool ok() const { return count > 0; }
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
 */
// ----------------------------------------------------------------------

void extrapolate(const std::unique_ptr<NFmiDataMatrix<float>> &theValues, int theAmount)
{
  if (theAmount <= 0)
    return;

  if (!theValues)
    throw Fmi::Exception(BCP, "Cannot extrapolate an empty matrix");

  auto &values = *theValues;

  const auto nan = std::numeric_limits<float>::quiet_NaN();

  auto tmp = values;
  for (int k = 0; k < theAmount; k++)
  {
    for (unsigned int j = 0; j < values.NY(); j++)
      for (unsigned int i = 0; i < values.NX(); i++)
      {
        if (std::isnan(tmp[i][j]))
        {
          Extrapolation ext;
          ext(values.At(i - 1, j, nan));
          ext(values.At(i + 1, j, nan));
          ext(values.At(i, j - 1, nan));
          ext(values.At(i, j + 1, nan));
          if (!ext.ok())
          {
            ext(values.At(i - 1, j - 1, nan));
            ext(values.At(i - 1, j + 1, nan));
            ext(values.At(i + 1, j - 1, nan));
            ext(values.At(i + 1, j + 1, nan));
          }
          if (ext.ok())
            tmp[i][j] = ext.result();
        }
      }

    std::swap(tmp, values);
  }
}

// ----------------------------------------------------------------------
/*!
 * \brief Flip values
 */
// ----------------------------------------------------------------------

void flip(NFmiDataMatrix<float> &values)
{
  const auto nx = values.NX();
  const auto ny = values.NY();

  for (std::size_t j1 = 0; j1 < ny / 2; j1++)
  {
    std::size_t j2 = ny - j1 - 1;

    for (std::size_t i = 0; i < nx; i++)
      std::swap(values[i][j1], values[i][j2]);
  }
}

// ----------------------------------------------------------------------
/*!
 * \brief Flip coordinates
 */
// ----------------------------------------------------------------------

void flip(Fmi::CoordinateMatrix &coords)
{
  const auto nx = coords.width();
  const auto ny = coords.height();

  for (std::size_t j1 = 0; j1 < ny / 2; j1++)
  {
    std::size_t j2 = ny - j1 - 1;

    for (std::size_t i = 0; i < nx; i++)
    {
      NFmiPoint pt1 = coords(i, j1);
      NFmiPoint pt2 = coords(i, j2);
      coords.set(i, j1, pt2);
      coords.set(i, j2, pt1);
    }
  }
}

// ----------------------------------------------------------------------
/*!
 * \brief Determine if a rectangular bounding box and a grid cell overlap
 *
 * We assume the bounding box to be much larger than grid cells, and hence
 * do not bother with special cases which can only be resolved by calculating
 * intersections. Knowing when an intersection is simply not possible is
 * enough for contouring purposes.
 */
// ----------------------------------------------------------------------

bool overlaps(double vmin, double vmax, double v1, double v2, double v3, double v4)
{
  // Everything to the right or above the bbox?
  if (v1 > vmax && v2 > vmax && v3 > vmax && v4 > vmax)
    return false;

  // Everything to the left or below the bbox?
  if (v1 < vmin && v2 < vmin && v3 < vmin && v4 < vmin)
    return false;

  // If a cell is marked as valid due to NaN coordinates, so be it. Other parts
  // of the code are expected to handle such cases.
  return true;
}

// ----------------------------------------------------------------------
/*!
 * \brief Mark cells which overlap the given bounding box
 */
// ----------------------------------------------------------------------

Fmi::BoolMatrix mark_valid_cells(const Fmi::CoordinateMatrix &theCoordinates,
                                 const Fmi::Box &theClipBox)
{
  // Assume no overlap for all cells
  const auto ny = theCoordinates.height();
  const auto nx = theCoordinates.width();

  Fmi::BoolMatrix ret(nx - 1, ny - 1, false);

  for (auto j = 0UL; j < ny - 1; j++)
    for (auto i = 0UL; i < nx - 1; i++)
    {
      if (overlaps(theClipBox.xmin(),
                   theClipBox.xmax(),
                   theCoordinates.x(i, j),
                   theCoordinates.x(i, j + 1),
                   theCoordinates.x(i + 1, j + 1),
                   theCoordinates.x(i + 1, j)) &&
          overlaps(theClipBox.ymin(),
                   theClipBox.ymax(),
                   theCoordinates.y(i, j),
                   theCoordinates.y(i, j + 1),
                   theCoordinates.y(i + 1, j + 1),
                   theCoordinates.y(i + 1, j)))
      {
        ret.set(i, j, true);
      }
    }

  return ret;
}

// ----------------------------------------------------------------------
/*!
 * \brief What coordinates to shift and by how much
 */
// ----------------------------------------------------------------------

struct GlobeShift
{
  GlobeShift() = delete;
  GlobeShift(double w, double e, double s) : west(w), east(e), shift(s) {}
  double west;
  double east;
  double shift;
  std::size_t hashValue() const
  {
    auto hash = Fmi::hash_value(west);
    Fmi::hash_combine(hash, Fmi::hash_value(east));
    Fmi::hash_combine(hash, Fmi::hash_value(shift));
    return hash;
  }
};

// ----------------------------------------------------------------------
/*!
 * \brief Get data box
 */
// ----------------------------------------------------------------------

Fmi::BBox get_bbox(const Fmi::CoordinateMatrix &theCoordinates)
{
  try
  {
    if (theCoordinates.width() == 0 || theCoordinates.height() == 0)
      throw Fmi::Exception(BCP, "Trying to extract BBOX of an empty grid");

    Fmi::BBox bbox;
    bool first = true;
    for (auto j = 0UL; j <= theCoordinates.height(); j++)
      for (auto i = 0UL; i <= theCoordinates.width(); i++)
      {
        const auto x = theCoordinates.x(i, j);
        const auto y = theCoordinates.y(i, j);
        if (!std::isnan(x) && !std::isnan(y))
        {
          if (first)
          {
            bbox.west = x;
            bbox.east = x;
            bbox.north = y;
            bbox.south = y;
            first = false;
          }
          else
          {
            bbox.west = std::min(bbox.west, x);
            bbox.east = std::max(bbox.east, x);
            bbox.south = std::min(bbox.south, y);
            bbox.north = std::max(bbox.north, y);
          }
        }
      }
    return bbox;
  }
  catch (...)
  {
    throw Fmi::Exception::Trace(BCP, "Operation failed!");
  }
}

// ----------------------------------------------------------------------
/*!
 * \brief Determine if the coordinates should be shifted to get a Pacific view
 */
// ----------------------------------------------------------------------

boost::optional<GlobeShift> get_globe_shift(const Fmi::CoordinateMatrix &theCoordinates,
                                            const Fmi::SpatialReference &theOutputCRS,
                                            const boost::optional<Fmi::BBox> &theBBox)
{
  // If there is a user requested BBOX, we need to check whether it is within the BBOX of the CRS.
  // If not, we will shift the coordinates by the width of the CRS BBOX to get for example a
  // Pacific view of the data. However, this only really works for global CRSs. Hence we also need
  // to verify the CRS WGS84 BBOX is global (-180...180, latitude does not matter). We could also
  // limit the shift to WebMercator (3857), but this solution is more general.

  if (!theBBox)
    return {};

  const auto epsg = theOutputCRS.getEPSG();
  if (!epsg)
    return {};

  auto info = Fmi::EPSGInfo::getInfo(*epsg);
  if (!info)
    return {};

  const auto &bbox = info->bbox;

  bool is_global = (bbox.west == -180 && bbox.east == 180);
  if (!is_global)
    return {};

  // Now we need to check the overlap. We allow for some tolerance due to
  // possible rounding errors in projection calculations due to software
  // differences. 0.05 degrees ~ 5000 meters, no need to be exact here.
  // We do this check first because calculating the bbox of the data takes time.

  const double tolerance = (theOutputCRS.isGeographic() ? 0.05 : 5000.0);

  const auto &bounds = info->bounds;
  const auto bounds_width = bounds.east - bounds.west;
  const auto bounds_center = (bounds.west + bounds.east) / 2;

  if (theBBox->west >= bounds.west - tolerance && theBBox->east <= bounds.east + tolerance)
    return {};

  // If the requested bbox fully contains the data, there is no need to shift anything.
  // For example we may have zoomed out so that the 180th meridian is crossed, but the
  // data might cover only Scandinavia and is hence still fully covered.

  const auto data_bbox = get_bbox(theCoordinates);
  if (data_bbox.west >= theBBox->west && data_bbox.east <= theBBox->east)
    return {};

  // Finally we actually generate a shift if necessary

  if (theBBox->west < bounds.west - tolerance)
    return GlobeShift{bounds_center, bounds.east, -bounds_width};  // right side of globe to left

  if (theBBox->east > bounds.east + tolerance)
    return GlobeShift{bounds.west, bounds_center, bounds_width};  // left side of globe to right

  return {};
}

std::unique_ptr<Fmi::CoordinateMatrix> shift_globe(const Fmi::CoordinateMatrix &coords,
                                                   const GlobeShift &globe_shift)
{
  std::unique_ptr<Fmi::CoordinateMatrix> ret(new Fmi::CoordinateMatrix(coords));

  for (auto j = 0UL; j < coords.height(); j++)
    for (auto i = 0UL; i < coords.width(); i++)
    {
      const auto x = coords.x(i, j);
      if (x >= globe_shift.west && x < globe_shift.east)
        ret->set(i, j, x + globe_shift.shift, coords.y(i, j));  // x,y --> x+width,y
    }
  return ret;
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
      throw Fmi::Exception(BCP, "Ill-defined isocircle: start and end points are equal");

    if (std::abs(lon1 - lon2) == 180 && std::abs(lat1 - (90 - lat2)) == 90)
      throw Fmi::Exception(BCP, "Ill-defined isocircle: points at opposing ends of the earth");

    if (steps < 1 || steps > 10000)
      throw Fmi::Exception(
          BCP, "Number of points on isocircle must be 1-10000, not " + Fmi::to_string(steps));

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
      const bool pacific_view = false;
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
    throw Fmi::Exception::Trace(BCP, "Operation failed!");
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
      return {};

    // Convert to WKB

    std::ostringstream out;
    geos::io::WKBWriter writer;
    writer.write(*theGeom, out);
    const auto &wkb = out.str();

    // Read to OGR using as hideous casts as is required

    auto *cwkb = reinterpret_cast<unsigned char *>(const_cast<char *>(wkb.c_str()));

    OGRGeometry *ogeom;
    OGRErr err = OGRGeometryFactory::createFromWkb(cwkb, theSR, &ogeom);

    if (err != OGRERR_NONE)
    {
      // we assume createFromWkb does not destroy the spatial reference on failure
      delete theSR;
      throw Fmi::Exception(BCP, "Failed to convert contoured WKB to OGRGeometry");
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
    throw Fmi::Exception::Trace(BCP, "Operation failed!");
  }
}

// ----------------------------------------------------------------------
/*!
 * \brief Implementation constructor
 */
// ----------------------------------------------------------------------

Engine::Impl::Impl(std::string theFileName)
    : itsConfigFile(std::move(theFileName)), itsGeomFactory(geos::geom::GeometryFactory::create())
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

    // Rough estimate: 50 producers times 10 projections = 500
    itsAnalysisCache.resize(1000);
  }
  catch (...)
  {
    throw Fmi::Exception::Trace(BCP, "Operation failed!");
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
    throw Fmi::Exception::Trace(BCP, "Operation failed!");
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
 * \brief Low level contouring utility for cross sections
 *
 * Note: We intentionally copy the values to be able to filter the data.
 *       We cache the data values, hence modifying the input is not OK.
 */
// ----------------------------------------------------------------------

GeometryPtr Engine::Impl::internal_isoline(const Trax::Grid &data,
                                           double isovalue,
                                           Trax::InterpolationType interpolation) const
{
  Trax::IsolineValues limits;
  limits.add(isovalue);

  Trax::Contour contourer;
  contourer.interpolation(interpolation);
  contourer.strict(false);
  contourer.validate(false);

  auto result = contourer.isolines(data, limits);
  return Trax::to_geos_geom(result[0], itsGeomFactory);
}

// ----------------------------------------------------------------------
/*!
 * \brief Low level contouring utility for cross sections
 *
 * Note: We intentionally copy the values to be able to filter the data.
 *       We cache the data values, hence modifying the input is not OK.
 */
// ----------------------------------------------------------------------

GeometryPtr Engine::Impl::internal_isoband(const Trax::Grid &data,
                                           const boost::optional<double> &lolimit,
                                           const boost::optional<double> &hilimit,
                                           Trax::InterpolationType interpolation) const
{
  // Change optional settings to corresponding contours. If either limit is set,
  // the other is +-inf if not set. If neither limit is set, we're contouring
  // missing values which for Trax is NaN...NaN. Thus if the user really wants
  // to contour all valid values, input -Inf...Inf must be given instead of
  // optional values.

  double lo = -std::numeric_limits<double>::infinity();
  double hi = +std::numeric_limits<double>::infinity();

  if (!lolimit && !hilimit)
  {
    lo = std::numeric_limits<double>::quiet_NaN();
    hi = lo;
  }
  else
  {
    if (lolimit)
      lo = *lolimit;
    if (hilimit)
      hi = *hilimit;
  }

  Trax::IsobandLimits limits;
  limits.add(lo, hi);

  Trax::Contour contourer;
  contourer.interpolation(interpolation);
  contourer.closed_range(true);
  contourer.strict(false);
  contourer.validate(false);
  auto result = contourer.isobands(data, limits);

  return Trax::to_geos_geom(result[0], itsGeomFactory);
}

// ----------------------------------------------------------------------
/*!
 * \brief Return a coordinate analysis of the projected grid coordinates
 */
// ----------------------------------------------------------------------

std::shared_ptr<Fmi::CoordinateAnalysis> Engine::Impl::get_analysis(
    std::size_t theDataHash,
    const Fmi::CoordinateMatrix &theCoordinates,
    const Fmi::SpatialReference &theOutputCRS) const
{
  // Combined hash value includes data details and output projection
  std::size_t hash = theDataHash;
  Fmi::hash_combine(hash, theOutputCRS.hashValue());

  // See if the result has already been cached
  auto cached = itsAnalysisCache.find(hash);
  if (cached)
    return cached->get();

  // Calculate the analysis asynchronously
  auto ftr =
      std::async(
          [&] { return std::make_shared<Fmi::CoordinateAnalysis>(Fmi::analysis(theCoordinates)); })
          .share();

  // Cache the future analysis to block other simlultaneous requests
  itsAnalysisCache.insert(hash, ftr);

  // And return the result when the async analysis finishes
  return ftr.get();
}

// ----------------------------------------------------------------------
/*!
 * \brief Calculate cache keys for the contours to be calculated
 */
// ----------------------------------------------------------------------

std::vector<std::size_t> get_contour_cache_keys(std::size_t theDataHash,
                                                const Fmi::SpatialReference &theOutputCRS,
                                                const Options &theOptions)
{
  // The hash for the result
  auto common_hash = theDataHash;
  Fmi::hash_combine(common_hash, theOptions.filtered_data_hash_value());
  Fmi::hash_combine(common_hash, theOutputCRS.hashValue());

  if (!theOptions.isovalues.empty())
  {
    // isolines mode
    const auto nresults = theOptions.isovalues.size();
    std::vector<std::size_t> hash_values(nresults, 0);

    const auto isoline_hash = Fmi::hash_value(std::string("isoline"));
    Fmi::hash_combine(common_hash, isoline_hash);

    for (auto i = 0UL; i < nresults; i++)
    {
      auto hash = common_hash;
      Fmi::hash_combine(hash, Fmi::hash_value(theOptions.isovalues[i]));
      hash_values[i] = hash;
    }

    return hash_values;
  }

  // isoband mode
  const auto nresults = theOptions.limits.size();
  std::vector<std::size_t> hash_values(nresults, 0);

  const auto isoband_hash = Fmi::hash_value(std::string("isoband"));
  Fmi::hash_combine(common_hash, isoband_hash);

  for (auto i = 0UL; i < nresults; i++)
  {
    auto hash = common_hash;

    // TODO: Use generic code as in wms/Hash.h
    if (theOptions.limits[i].lolimit)
      Fmi::hash_combine(hash, Fmi::hash_value(*theOptions.limits[i].lolimit));
    else
      Fmi::hash_combine(hash, Fmi::hash_value(false));

    if (theOptions.limits[i].hilimit)
      Fmi::hash_combine(hash, Fmi::hash_value(*theOptions.limits[i].hilimit));
    else
      Fmi::hash_combine(hash, Fmi::hash_value(false));

    hash_values[i] = hash;
  }

  return hash_values;
}

// ----------------------------------------------------------------------
/*!
 * \brief Contour producing vector of OGR geometries
 */
// ----------------------------------------------------------------------

std::vector<OGRGeometryPtr> Engine::Impl::contour(std::size_t theDataHash,
                                                  const Fmi::SpatialReference &theOutputCRS,
                                                  const NFmiDataMatrix<float> &theMatrix,
                                                  const Fmi::CoordinateMatrix &theCoordinates,
                                                  const Options &theOptions) const
{
#if 0  
  // All valid cells are to be contoured
  Fmi::BoolMatrix valid_cells(theCoordinates.width() - 1, theCoordinates.height() - 1, true);

  // Dummy clipbox for generating a rectangle around the valid areas when contouring missing values
  const double large = 1E10;
  Fmi::Box clipbox(-large, -large, large, large, 1, 1);

  return contour(theDataHash, theOutputCRS, theMatrix, theCoordinates, valid_cells, theOptions);
#else
  // Dummy clipbox for generating a rectangle around the valid areas when contouring missing values
  const double large = 1E10;
  Fmi::Box clipbox(-large, -large, large, large, 1, 1);

  return contour(theDataHash, theOutputCRS, theMatrix, theCoordinates, clipbox, false, theOptions);
#endif
}

std::vector<OGRGeometryPtr> Engine::Impl::contour(std::size_t theDataHash,
                                                  const Fmi::SpatialReference &theOutputCRS,
                                                  const NFmiDataMatrix<float> &theMatrix,
                                                  const Fmi::CoordinateMatrix &theCoordinates,
                                                  const Fmi::Box &theClipBox,
                                                  const Options &theOptions) const
{
#if 0
  // Mark cells overlapping the bounding box
  auto valid_cells = mark_valid_cells(theCoordinates, theClipBox);
  return contour(theDataHash, theOutputCRS, theMatrix, theCoordinates, valid_cells, theOptions);
#else
  return contour(
      theDataHash, theOutputCRS, theMatrix, theCoordinates, theClipBox, true, theOptions);
#endif
}

std::vector<OGRGeometryPtr> Engine::Impl::contour(std::size_t theDataHash,
                                                  const Fmi::SpatialReference &theOutputCRS,
                                                  const NFmiDataMatrix<float> &theMatrix,
                                                  const Fmi::CoordinateMatrix &theCoordinates,
                                                  const Fmi::Box &theClipBox,
                                                  bool all_valid,
                                                  Options theOptions) const
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

    if (!theOptions.isovalues.empty() && !theOptions.limits.empty())
      throw Fmi::Exception(BCP, "Cannot calculate isolines and isobands simultaneously");

    // Change the contour limits with reverse unit conversion if required. This must be done
    // before the contour cache keys are calculated below, or the unit conversion has no
    // effect on the cache key.

    theOptions.transform();

    // Calculate the cache keys for all the contours. This enables us to test
    // if everything is cached already without doing any data processing.

    auto contour_hash = theDataHash;
    Fmi::hash_combine(contour_hash, theClipBox.hashValue());

    auto contour_cache_keys = get_contour_cache_keys(contour_hash, theOutputCRS, theOptions);

    // Initialize empty results

    auto nresults = std::max(theOptions.isovalues.size(), theOptions.limits.size());
    std::vector<OGRGeometryPtr> retval(nresults);

    // Try to extract everything from the cache first without any
    // preprocessing or parallelism

    std::vector<std::size_t> todo_contours;

    for (auto i = 0UL; i < nresults; i++)
    {
      auto hash = contour_cache_keys[i];
      auto opt_geom = itsContourCache.find(hash);
      if (opt_geom)
        retval[i] = *opt_geom;
      else
        todo_contours.push_back(i);
    }

    // We're done if there is nothing on the TODO-list

    if (todo_contours.empty())
      return retval;

    // See if the global data needs to be shifted

    const auto *coordinates = &theCoordinates;               // assuming not
    std::unique_ptr<Fmi::CoordinateMatrix> alt_coordinates;  // holder for shifted coordinates
    auto coordinates_hash = theCoordinates.hashValue();      // hash value for coordinates

    auto globe_shift = get_globe_shift(theCoordinates, theOutputCRS, theOptions.bbox);

    if (globe_shift)
    {
      alt_coordinates = shift_globe(theCoordinates, *globe_shift);    // shifted coordinates
      coordinates = alt_coordinates.get();                            // new active coordinates
      Fmi::hash_combine(coordinates_hash, globe_shift->hashValue());  // new active hash
    }
    // From now on we use coordinates* instead of theCoordinates

    // Get the grid analysis of the projected coordinates from
    // the cache, or do the analysis and cache it.

    auto analysis = get_analysis(coordinates_hash, *coordinates, theOutputCRS);

    /*
     * Finally we determine which grid cells are valid and have the correct winding rule
     * while taking into account potential flipping.
     *
     *    Needs flipping:
     *         F   T
     *       +-------
     * CW  T | T   F      --> xor returns correct combination of CW cells even when flipped
     *     F | F   T
     *
     * Note: Fmi::analysis could do this for us, and the result would be cached. This
     *       is very fast though.
     *
     * Note: Contouring may be tiled in which case 'theValidCells' marks which cell overlap
     *       with the bouding box (expanded by clipping margin). We combine the results
     *       here instead of optimizing the analysis by skipping cells outside the tile
     *       since the analysis is cached for the full grid and hence there are no real
     *       speed gains to be obtained for WMS services.
     */

    auto valid_cells = analysis->valid;
    valid_cells &= (analysis->clockwise ^ analysis->needs_flipping);
    if (!all_valid)
      valid_cells &= mark_valid_cells(*coordinates, theClipBox);

    // Process the data for contouring. We wish to avoid unnecessary copying of
    // the data, hence we use the existence of an alternative unique_ptr
    // to indicate the data has been processed.

    std::unique_ptr<NFmiDataMatrix<float>> alt_values;

    const bool needs_copying = (theOptions.extrapolation > 0 || analysis->needs_flipping ||
                                theOptions.filter_size || theOptions.filter_degree);

    if (needs_copying)
      alt_values.reset(new NFmiDataMatrix<float>(theMatrix));

    // Extrapolate if requested
    extrapolate(alt_values, theOptions.extrapolation);

    // Flip the data and the coordinates if necessary. We avoid an unnecessary
    // copy of the coordinates if no flipping is needed by using a pointer.

    if (analysis->needs_flipping)
    {
      flip(*alt_values);
      alt_coordinates.reset(new Fmi::CoordinateMatrix(*coordinates));
      flip(*alt_coordinates);
      coordinates = alt_coordinates.get();
    }

    // Helper data structure for contouring. To be updated to std::unique_ptr with C++17
    std::shared_ptr<BaseGrid> data;

    if (analysis->shift == 0)
    {
      if (alt_values)
        data = std::make_shared<PaddedGrid>(*alt_values, *coordinates, valid_cells);
      else
      {
        // We're not really modifying the data, only alt_values is
        auto &tmp = const_cast<NFmiDataMatrix<float> &>(theMatrix);
        data = std::make_shared<PaddedGrid>(tmp, *coordinates, valid_cells);
      }
    }
    else
    {
      if (alt_values)
        data =
            std::make_shared<ShiftedGrid>(*alt_values, *coordinates, valid_cells, analysis->shift);
      else
      {
        // We're not really modifying the data, only alt_values is
        auto &tmp = const_cast<NFmiDataMatrix<float> &>(theMatrix);
        data = std::make_shared<ShiftedGrid>(tmp, *coordinates, valid_cells, analysis->shift);
      }
    }

    // Savitzky-Golay assumes DataMatrix like Adapter API, not NFmiDataMatrix like API, hence the
    // filtering is delayed until this point.

    if (theOptions.filter_size || theOptions.filter_degree)
    {
      size_t size = (theOptions.filter_size ? *theOptions.filter_size : 1);
      size_t degree = (theOptions.filter_degree ? *theOptions.filter_degree : 1);
      data->smooth(size, degree);
    }

    Trax::Contour contourer;
    contourer.interpolation(theOptions.interpolation);
    contourer.closed_range(theOptions.closed_range);
    contourer.strict(theOptions.strict);
    contourer.validate(theOptions.validate);
    contourer.desliver(theOptions.desliver);

    Trax::GeometryCollections results;

    if (!theOptions.isovalues.empty())
    {
      // Calculate isolines
      Trax::IsolineValues isovalues;
      for (auto i : todo_contours)
        isovalues.add(theOptions.isovalues[i]);

      results = contourer.isolines(*data, isovalues);
    }
    else
    {
      // Calculate isobands
      Trax::IsobandLimits isobands;
      for (auto i : todo_contours)
      {
        // Change optional settings to corresponding contours. If either limit is set,
        // the other is +-inf if not set. If neither limit is set, we're contouring
        // missing values which for Trax is NaN...NaN. Thus if the user really wants
        // to contour all valid values, input -Inf...Inf must be given instead of
        // optional values.
        double lo = -std::numeric_limits<double>::infinity();
        double hi = +std::numeric_limits<double>::infinity();
        const auto &limits = theOptions.limits[i];

        if (!limits.lolimit && !limits.hilimit)
        {
          lo = std::numeric_limits<double>::quiet_NaN();
          hi = lo;
        }
        else
        {
          if (limits.lolimit)
            lo = *limits.lolimit;
          if (limits.hilimit)
            hi = *limits.hilimit;
        }

        isobands.add(lo, hi);
      }

      results = contourer.isobands(*data, isobands);
    }

    // Update results and cache
    for (auto i = 0UL; i < todo_contours.size(); i++)
    {
      OGRGeometryPtr tmp(Trax::to_ogr_geom(results[i]).release());

      if (theOutputCRS.get() != nullptr)
	tmp->assignSpatialReference(theOutputCRS.get()); // increments ref count

      // Despeckle even closed isolines (pressure curves)
      if (theOptions.minarea)
        tmp.reset(Fmi::OGR::despeckle(*tmp, *theOptions.minarea));

      retval[todo_contours[i]] = tmp;
      itsContourCache.insert(contour_cache_keys[todo_contours[i]], tmp);
    }

    return retval;
  }
  catch (...)
  {
    throw Fmi::Exception::Trace(BCP, "Operation failed!");
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
      throw Fmi::Exception(BCP, "Using the level option is meaningless for cross-sections");
    if (theOptions.filter_size)
      throw Fmi::Exception(BCP, "Using the filter_size option is meaningless for cross-sections");
    if (theOptions.filter_degree)
      throw Fmi::Exception(BCP, "Using the filter_degree option is meaningless for cross-sections");

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

        auto value = theQInfo.InterpolatedValue(coordinates[i], theOptions.time, max_minutes);
        if (value == kFloatMissing)
          value = std::numeric_limits<float>::quiet_NaN();

        values[i][j] = value;

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

    // Don't bother analyzing the grid for cross sections, the Z coordinates
    // should always be fine.

    Fmi::BoolMatrix grid_is_fine(coords.width() - 1, coords.height() - 1, true);
    PaddedGrid data(values, coords, grid_is_fine);

    // results first include isolines, then isobands
    auto nresults = theOptions.isovalues.size() + theOptions.limits.size();

    // The results
    std::vector<OGRGeometryPtr> retval(nresults);

    for (auto icontour = 0UL; icontour < nresults; icontour++)
    {
      // Is it an isovalue request?
      bool isovaluerequest = (icontour < theOptions.isovalues.size());

      // Calculate GEOS geometry result
      GeometryPtr geom;
      if (isovaluerequest)
      {
        geom = internal_isoline(data, theOptions.isovalues[icontour], theOptions.interpolation);
      }
      else
      {
        auto i = icontour - theOptions.isovalues.size();
        geom = internal_isoband(data,
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
    throw Fmi::Exception::Trace(BCP, "Operation failed!");
  }
}

Fmi::Cache::CacheStatistics Engine::Impl::getCacheStats() const
{
  Fmi::Cache::CacheStatistics ret;

  ret.insert(std::make_pair("Contour::contour_cache", itsContourCache.statistics()));
  ret.insert(std::make_pair("Contour::analysis_cache", itsAnalysisCache.statistics()));

  return ret;
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
    throw Fmi::Exception::Trace(BCP, "Operation failed!");
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
    throw Fmi::Exception::Trace(BCP, "Operation failed!");
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

std::vector<OGRGeometryPtr> Engine::contour(std::size_t theDataHash,
                                            const Fmi::SpatialReference &theOutputCRS,
                                            const NFmiDataMatrix<float> &theMatrix,
                                            const Fmi::CoordinateMatrix &theCoordinates,
                                            const Options &theOptions) const
{
  try
  {
    return itsImpl->contour(theDataHash, theOutputCRS, theMatrix, theCoordinates, theOptions);
  }
  catch (...)
  {
    throw Fmi::Exception::Trace(BCP, "Operation failed!");
  }
}

// ----------------------------------------------------------------------
/*!
 * \brief Contour given Box only
 */
// ----------------------------------------------------------------------

std::vector<OGRGeometryPtr> Engine::contour(std::size_t theDataHash,
                                            const Fmi::SpatialReference &theOutputCRS,
                                            const NFmiDataMatrix<float> &theMatrix,
                                            const Fmi::CoordinateMatrix &theCoordinates,
                                            const Fmi::Box &theClipBox,
                                            const Options &theOptions) const
{
  try
  {
    return itsImpl->contour(
        theDataHash, theOutputCRS, theMatrix, theCoordinates, theClipBox, theOptions);
  }
  catch (...)
  {
    throw Fmi::Exception::Trace(BCP, "Operation failed!");
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
    throw Fmi::Exception::Trace(BCP, "Operation failed!");
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
    throw Fmi::Exception::Trace(BCP, "Operation failed!");
  }
}

Fmi::Cache::CacheStatistics Engine::getCacheStats() const
{
  return itsImpl->getCacheStats();
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
