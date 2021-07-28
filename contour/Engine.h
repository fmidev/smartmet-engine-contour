// ======================================================================
/*!
 * \brief Interface of class Contour
 */
// ======================================================================

#pragma once

#include "Options.h"
#include <gis/CoordinateMatrix.h>
#include <gis/SpatialReference.h>
#include <gis/Types.h>
#include <newbase/NFmiFastQueryInfo.h>
#include <newbase/NFmiParameterName.h>
#include <spine/SmartMetEngine.h>
#include <memory>
#include <string>

class OGRGeometry;
class OGRSpatialReference;

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
 public:
  // constructor is available only with a libconfig configuration file

  Engine() = delete;
  Engine(const std::string& theFileName);

  // Produce vector of OGR geometries in output spatial reference

  std::vector<OGRGeometryPtr> contour(std::size_t theDataHash,
                                      const Fmi::SpatialReference& theDataCRS,
                                      const Fmi::SpatialReference& theOutputCRS,
                                      const NFmiDataMatrix<float>& theMatrix,
                                      const Fmi::CoordinateMatrix& theCoordinates,
                                      const Options& theOptions) const;

  // Produce a cross section contour
  std::vector<OGRGeometryPtr> crossection(NFmiFastQueryInfo& theQInfo,
                                          const Options& theOptions,
                                          double theLon1,
                                          double theLat1,
                                          double theLon2,
                                          double theLat2,
                                          std::size_t theSteps) const;

  // Produce a cross section contour with given parameter for Z-values
  std::vector<OGRGeometryPtr> crossection(NFmiFastQueryInfo& theQInfo,
                                          const Spine::Parameter& theZParameter,
                                          const Options& theOptions,
                                          double theLon1,
                                          double theLat1,
                                          double theLon2,
                                          double theLat2,
                                          std::size_t theSteps) const;

  CacheReportingStruct getCacheSizes() const;
  void clearCache() const;

 protected:
  void init() override;
  void shutdown() override;

 private:
  // Hide caches etc behind an implementation
  class Impl;
  std::shared_ptr<Impl> itsImpl;

};  // class Engine

}  // namespace Contour
}  // namespace Engine
}  // namespace SmartMet

// ======================================================================
