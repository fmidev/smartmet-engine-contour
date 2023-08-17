#pragma once

#include "NormalGrid.h"

namespace SmartMet
{
namespace Engine
{
namespace Contour
{
namespace SavitzkyGolay2D
{
// ----------------------------------------------------------------------
/*!
 * \brief Smoothen a matrix using a mirror matrix for boundary conditions
 *
 * Size is limited to range 0..6, degree to 0...5
 */
// ----------------------------------------------------------------------

void smooth(NormalGrid& input, std::size_t size, std::size_t degree);

}  // namespace SavitzkyGolay2D
}  // namespace Contour
}  // namespace Engine
}  // namespace SmartMet
