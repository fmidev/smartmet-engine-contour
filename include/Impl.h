// ======================================================================
/*!
 * \brief Implementation details for the Contour engine
 *
 * TODO: Should have a separate GeometryFactory for each unique SRID?
 */
// ======================================================================

#pragma once

#include "Config.h"
#include "DataMatrixAdapter.h"
#include "Engine.h"

#include <tron/FmiBuilder.h>
#include <tron/Tron.h>
#include <macgyver/Cache.h>
#include <geos/geom/Geometry.h>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

typedef boost::shared_ptr<geos::geom::Geometry> GeometryPtr;

// Contourers

typedef Tron::Traits<double, double, Tron::FmiMissing> MyTraits;

typedef Tron::Contourer<DataMatrixAdapter, Tron::FmiBuilder, MyTraits, Tron::LinearInterpolation>
    MyLinearContourer;

typedef Tron::Contourer<DataMatrixAdapter, Tron::FmiBuilder, MyTraits, Tron::LogLinearInterpolation>
    MyLogLinearContourer;

typedef Tron::Contourer<DataMatrixAdapter,
                        Tron::FmiBuilder,
                        MyTraits,
                        Tron::NearestNeighbourInterpolation> MyNearestContourer;

typedef Tron::Contourer<DataMatrixAdapter, Tron::FmiBuilder, MyTraits, Tron::DiscreteInterpolation>
    MyDiscreteContourer;

typedef Tron::Hints<DataMatrixAdapter, MyTraits> MyHints;

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
  void init();

  // Produce vector of OGR contours for the given spatial reference
  std::vector<OGRGeometryPtr> contour(std::size_t theQhash,
                                      const std::string theQAreaWKT,
                                      const NFmiDataMatrix<float> &theMatrix,
                                      const CoordinatesPtr theCoordinates,
                                      const Options &theOptions,
                                      OGRSpatialReference *theSR);

  // Produce an OGR crossection for the given data
  std::vector<OGRGeometryPtr> crossection(
      boost::shared_ptr<NFmiFastQueryInfo> theQInfo,
      const Options &theOptions,
      const boost::optional<SmartMet::Spine::Parameter> &theZParameter,
      double theLon1,
      double theLat1,
      double theLon2,
      double theLat2,
      std::size_t theSteps);

  CacheReportingStruct getCacheSizes();

 private:
  // Private parts

  Impl();
  std::string itsConfigFile;
  std::unique_ptr<Config> itsConfig;  // ptr for delayed initialization

  boost::shared_ptr<geos::geom::GeometryFactory> itsGeomFactory;

  GeometryPtr internal_contour(std::size_t datahash,
                               const Options &options,
                               const DataMatrixAdapter &data,
                               const MyHints &hints);
  // Cached contours

  typedef Fmi::Cache::Cache<std::size_t, OGRGeometryPtr> GeometryCache;
  mutable GeometryCache itsContourCache;

  GeometryPtr geosContour(std::size_t theQhash,
                          const Options &theOptions,
                          const DataMatrixAdapter &theData,
                          const MyHints &theHints);
};

}  // namespace Contour
}  // namespace Engine
}  // namespace SmartMet
