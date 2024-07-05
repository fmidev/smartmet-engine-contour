#include "Config.h"
#include <filesystem>
#include <macgyver/Exception.h>
#include <spine/ConfigTools.h>
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

Config::Config(const std::string& theFileName)
{
  try
  {
    if (theFileName.empty())
      throw Fmi::Exception(BCP, "Contour-engine config filename is empty");

    // Enable sensible relative include paths
    std::filesystem::path p = theFileName;
    p.remove_filename();
    itsConfig.setIncludeDir(p.c_str());

    itsConfig.readFile(theFileName.c_str());
    Spine::expandVariables(itsConfig);

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
