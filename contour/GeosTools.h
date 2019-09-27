// ======================================================================
/*!
 * \brief Various GEOS related tools
 */
// ======================================================================

#pragma once

#include <string>

namespace geos
{
namespace geom
{
class Geometry;
}
}  // namespace geos

namespace SmartMet
{
namespace GeosTools
{
std::string getSVG(const geos::geom::Geometry& theGeom, int thePrecision = -1);
}
}  // namespace SmartMet
