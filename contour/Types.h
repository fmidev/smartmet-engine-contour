// ======================================================================
/*!
 * \brief Common type definitions
 */
// ======================================================================

#pragma once
#include <gdal/ogr_geometry.h>
#include <memory>

namespace SmartMet
{
typedef std::shared_ptr<OGRGeometry> OGRGeometryPtr;
}
