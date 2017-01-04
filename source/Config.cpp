#include "Config.h"
#include <spine/Exception.h>
#include <stdexcept>

namespace SmartMet
{
namespace Engine
{
namespace Contour
{
// ----------------------------------------------------------------------
/*!
 * \brief Config constructor
 */
// ----------------------------------------------------------------------

Config::Config(const std::string& theFileName) : itsConfig()
{
  try
  {
    if (theFileName.empty())
      throw SmartMet::Spine::Exception(BCP, "Contour-engine config filename is empty");

    itsConfig.readFile(theFileName.c_str());

    if (!itsConfig.exists("cache.max_contours"))
      throw SmartMet::Spine::Exception(BCP, "cache.max_contours not set in '" + theFileName + "'");
    itsConfig.lookupValue("cache.max_contours", itsMaxContourCacheSize);
  }
  catch (...)
  {
    throw SmartMet::Spine::Exception(BCP, "Operation failed!", NULL);
  }
}

}  // namespace Contour
}  // namespace Engine
}  // namespace SmartMet