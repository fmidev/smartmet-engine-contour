// ======================================================================
/*!
 * \brief Common type definitions
 */
// ======================================================================

#pragma once
#include <boost/shared_ptr.hpp>
#include <ogr_geometry.h>

namespace SmartMet
{
typedef boost::shared_ptr<OGRGeometry> OGRGeometryPtr;
}
