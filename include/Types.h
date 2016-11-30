// ======================================================================
/*!
 * \brief Common type definitions
 */
// ======================================================================

#pragma once
#include <gdal/ogr_geometry.h>
#include <boost/shared_ptr.hpp>

namespace SmartMet
{
typedef boost::shared_ptr<OGRGeometry> OGRGeometryPtr;
}
