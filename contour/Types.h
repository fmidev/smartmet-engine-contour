// ======================================================================
/*!
 * \brief Common type definitions
 */
// ======================================================================

#pragma once
#include <gdal/ogr_geometry.h>
#include <shared_ptr>

namespace SmartMet
{
typedef std::shared_ptr<OGRGeometry> OGRGeometryPtr;
}
