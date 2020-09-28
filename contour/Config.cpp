#include "Config.h"
#include <macgyver/Exception.h>
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
      throw Fmi::Exception(BCP, "Contour-engine config filename is empty");

    itsConfig.readFile(theFileName.c_str());

    if (!itsConfig.exists("cache.max_contours"))
      throw Fmi::Exception(BCP, "cache.max_contours not set in '" + theFileName + "'");
    itsConfig.lookupValue("cache.max_contours", itsMaxContourCacheSize);
  }
  catch (...)
  {
    throw Fmi::Exception::Trace(BCP, "Operation failed!");
  }
}

}  // namespace Contour
}  // namespace Engine
}  // namespace SmartMet
