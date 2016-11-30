// ======================================================================
/*!
 * \brief Implementation of Contour::Engine
 */
// ======================================================================

#include "Engine.h"
#include "Impl.h"
#include "GeosTools.h"
#include <spine/Exception.h>

#include <gdal/cpl_conv.h>  // For configuring GDAL

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

Engine::Engine(const std::string& theFileName) : itsImpl(new Impl(theFileName))
{
}
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
    throw SmartMet::Spine::Exception(BCP, "Operation failed!", NULL);
  }
}

CacheReportingStruct Engine::getCacheSizes() const
{
  try
  {
    return itsImpl->getCacheSizes();
  }
  catch (...)
  {
    throw SmartMet::Spine::Exception(BCP, "Operation failed!", NULL);
  }
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
                                            const std::string theQAreaWKT,
                                            const NFmiDataMatrix<float>& theMatrix,
                                            const CoordinatesPtr theCoordinates,
                                            const Options& theOptions,
                                            OGRSpatialReference* theSR) const
{
  try
  {
    return itsImpl->contour(theQhash, theQAreaWKT, theMatrix, theCoordinates, theOptions, theSR);
  }
  catch (...)
  {
    throw SmartMet::Spine::Exception(BCP, "Operation failed!", NULL);
  }
}

// ----------------------------------------------------------------------
/*!
 * \brief Cross section
 */
// ----------------------------------------------------------------------

std::vector<OGRGeometryPtr> Engine::crossection(boost::shared_ptr<NFmiFastQueryInfo> theQInfo,
                                                const Options& theOptions,
                                                double theLon1,
                                                double theLat1,
                                                double theLon2,
                                                double theLat2,
                                                std::size_t theSteps) const
{
  try
  {
    boost::optional<SmartMet::Spine::Parameter> zparam;
    return itsImpl->crossection(
        theQInfo, theOptions, zparam, theLon1, theLat1, theLon2, theLat2, theSteps);
  }
  catch (...)
  {
    throw SmartMet::Spine::Exception(BCP, "Operation failed!", NULL);
  }
}

// ----------------------------------------------------------------------
/*!
 * \brief Cross section
 */
// ----------------------------------------------------------------------

std::vector<OGRGeometryPtr> Engine::crossection(boost::shared_ptr<NFmiFastQueryInfo> theQInfo,
                                                const SmartMet::Spine::Parameter& theZParameter,
                                                const Options& theOptions,
                                                double theLon1,
                                                double theLat1,
                                                double theLon2,
                                                double theLat2,
                                                std::size_t theSteps) const
{
  try
  {
    boost::optional<SmartMet::Spine::Parameter> zparam(theZParameter);
    return itsImpl->crossection(
        theQInfo, theOptions, zparam, theLon1, theLat1, theLon2, theLat2, theSteps);
  }
  catch (...)
  {
    throw SmartMet::Spine::Exception(BCP, "Operation failed!", NULL);
  }
}

}  // namespace Contour
}  // namespace Engine
}  // namespace SmartMet

// DYNAMIC MODULE CREATION TOOLS

extern "C" void* engine_class_creator(const char* configfile, void* /* user_data */)
{
  return new SmartMet::Engine::Contour::Engine(configfile);
}

extern "C" const char* engine_name()
{
  return "Contour";
}
// ======================================================================
