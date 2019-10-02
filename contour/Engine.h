// ======================================================================
/*!
 * \brief Interface of class Contour
 */
// ======================================================================

#pragma once

#include "Options.h"
#include "Types.h"
#include <gdal/ogr_geometry.h>
#include <gdal/ogr_spatialref.h>
#include <newbase/NFmiFastQueryInfo.h>
#include <newbase/NFmiParameterName.h>
#include <spine/SmartMetEngine.h>
#include <shared_ptr>
#include <string>

using Coordinates = NFmiDataMatrix<NFmiPoint>;
using CoordinatesPtr = std::shared_ptr<Coordinates>;

namespace SmartMet
{
namespace Engine
{
namespace Contour
{
struct CacheReportingStruct
{
  int contour_cache_max_size;
  int contour_cache_size;
};

class Engine : public Spine::SmartMetEngine
{
 private:
  Engine();

  // Hide caches etc behind an implementation
  class Impl;
  std::shared_ptr<Impl> itsImpl;

 protected:
  virtual void init();
  void shutdown();

 public:
  // constructor is available only with a libconfig configuration file

  Engine(const std::string& theFileName);

  // Produce vector of OGR geometries in world XY coordinates

  std::vector<OGRGeometryPtr> contour(std::size_t theQhash,
                                      const std::string theQAreaWKT,
                                      const NFmiDataMatrix<float>& theMatrix,
                                      const CoordinatesPtr theCoordinates,
                                      const Options& theOptions,
                                      bool worldwrap,
                                      OGRSpatialReference* theSR = 0) const;

  // Produce a cross section contour
  std::vector<OGRGeometryPtr> crossection(std::shared_ptr<NFmiFastQueryInfo> theQInfo,
                                          const Options& theOptions,
                                          double theLon1,
                                          double theLat1,
                                          double theLon2,
                                          double theLat2,
                                          std::size_t theSteps) const;

  // Produce a cross section contour with given parameter for Z-values
  std::vector<OGRGeometryPtr> crossection(std::shared_ptr<NFmiFastQueryInfo> theQInfo,
                                          const Spine::Parameter& theZParameter,
                                          const Options& theOptions,
                                          double theLon1,
                                          double theLat1,
                                          double theLon2,
                                          double theLat2,
                                          std::size_t theSteps) const;

  CacheReportingStruct getCacheSizes() const;
  void clearCache() const;

};  // class Engine

}  // namespace Contour
}  // namespace Engine
}  // namespace SmartMet

// ======================================================================
