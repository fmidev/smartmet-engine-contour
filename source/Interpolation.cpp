// ======================================================================
/*!
 * \brief Contour interpolation methods
 */
// ======================================================================

#include "Interpolation.h"
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
 * \brief Parse contour interpolation method name
 */
// ----------------------------------------------------------------------

Interpolation parseInterpolation(const std::string& theName)
{
  try
  {
    if (theName == "Nearest")
      return Nearest;
    else if (theName == "Linear")
      return Linear;
    else if (theName == "LogLinear")
      return LogLinear;
    else if (theName == "Discrete")
      return Discrete;

    throw SmartMet::Spine::Exception(BCP, "Unknown interpolation method '" + theName + "'");
  }
  catch (...)
  {
    throw SmartMet::Spine::Exception(BCP, "Operation failed!", NULL);
  }
}

}  // namespace Contour
}  // namespace Engine
}  // namespace SmartMet
