// ======================================================================
/*!
 * \brief Contour interpolation methods
 */
//======================================================================

#pragma once

#include <string>

namespace SmartMet
{
namespace Engine
{
namespace Contour
{
enum Interpolation
{
  Nearest,
  Linear,
  Discrete,
  LogLinear
};

Interpolation parseInterpolation(const std::string& theName);

}  // namespace Contour
}  // namespace Engine
}  // namespace SmartMet
