// ======================================================================
/*!
 * \brief Geos tools
 */
// ======================================================================

#include "GeosTools.h"
#include <boost/numeric/conversion/cast.hpp>
#include <geos/geom/Coordinate.h>
#include <geos/geom/CoordinateSequence.h>
#include <geos/geom/LineString.h>
#include <geos/geom/LinearRing.h>
#include <geos/geom/MultiLineString.h>
#include <geos/geom/MultiPoint.h>
#include <geos/geom/MultiPolygon.h>
#include <geos/geom/Point.h>
#include <geos/geom/Polygon.h>
#include <geos/geom/PrecisionModel.h>
#include <macgyver/Exception.h>
#include <iomanip>
#include <sstream>
#include <stdexcept>

using namespace geos::geom;

namespace SmartMet
{
namespace
{
void writeSVG(const Geometry* geom, std::ostream& out);

// ----------------------------------------------------------------------
/*!
 * \brief Handle a coordinate
 */
// ----------------------------------------------------------------------

void writeCoordinateSVG(const Coordinate& geom, std::ostream& out)
{
  out << geom.x << ' ' << geom.y;
}

// ----------------------------------------------------------------------
/*!
 * \brief Handle a Point
 */
// ----------------------------------------------------------------------

void writePointSVG(const Coordinate* geom, std::ostream& out)
{
  if (geom != nullptr)
    out << 'M' << geom->x << ' ' << geom->y;
}

// ----------------------------------------------------------------------
/*!
 * \brief Handle LinearRing
 */
// ----------------------------------------------------------------------

void writeLinearRingSVG(const LinearRing* geom, std::ostream& out)
{
  try
  {
    if (geom == nullptr || geom->isEmpty())
      return;
    // Bizarre: n is unsigned long but getting coordinate uses int
    for (unsigned long i = 0, n = geom->getNumPoints(); i < n - 1; ++i)
    {
      if (i == 0)
        out << 'M';
      else if (i == 1)
        out << ' ';  // implicit lineto, we could also write an 'L'
      else
        out << ' ';
      writeCoordinateSVG(geom->getCoordinateN(boost::numeric_cast<int>(i)), out);
    }
    out << 'Z';
  }
  catch (...)
  {
    throw Fmi::Exception::Trace(BCP, "Operation failed!");
  }
}

// ----------------------------------------------------------------------
/*!
 * \brief Handle LineString
 */
// ----------------------------------------------------------------------

void writeLineStringSVG(const LineString* geom, std::ostream& out)
{
  try
  {
    if (geom == nullptr || geom->isEmpty())
      return;

    unsigned long n = geom->getNumPoints();
    for (unsigned long i = 0; i < n - 1; ++i)
    {
      if (i == 0)
        out << 'M';
      else if (i == 1)
        out << ' ';  // implicit lineto, we could also write an 'L'
      else
        out << ' ';
      writeCoordinateSVG(geom->getCoordinateN(boost::numeric_cast<int>(i)), out);
    }

    if (geom->isClosed())
      out << 'Z';
    else
    {
      out << (n == 1 ? 'M' : ' ');
      writeCoordinateSVG(geom->getCoordinateN(boost::numeric_cast<int>(n - 1)), out);
    }
  }
  catch (...)
  {
    throw Fmi::Exception::Trace(BCP, "Operation failed!");
  }
}

// ----------------------------------------------------------------------
/*!
 * \brief Handle Polygon
 */
// ----------------------------------------------------------------------

void writePolygonSVG(const Polygon* geom, std::ostream& out)
{
  try
  {
    if (geom == nullptr || geom->isEmpty())
      return;

    writeLineStringSVG(geom->getExteriorRing(), out);
    for (size_t i = 0, n = geom->getNumInteriorRing(); i < n; ++i)
    {
      writeLineStringSVG(geom->getInteriorRingN(i), out);
    }
  }
  catch (...)
  {
    throw Fmi::Exception::Trace(BCP, "Operation failed!");
  }
}

// ----------------------------------------------------------------------
/*!
 * \brief Handle MultiPoint
 */
// ----------------------------------------------------------------------

void writeMultiPointSVG(const MultiPoint* geom, std::ostream& out)
{
  try
  {
    if (geom == nullptr || geom->isEmpty())
      return;

    for (size_t i = 0, n = geom->getNumGeometries(); i < n; ++i)
    {
      out << 'M';
      writeCoordinateSVG(*geom->getGeometryN(i)->getCoordinate(), out);
    }
  }
  catch (...)
  {
    throw Fmi::Exception::Trace(BCP, "Operation failed!");
  }
}

// ----------------------------------------------------------------------
/*!
 * \brief Handle MultiLineString
 */
// ----------------------------------------------------------------------

void writeMultiLineStringSVG(const MultiLineString* geom, std::ostream& out)
{
  try
  {
    if (geom == nullptr || geom->isEmpty())
      return;
    for (size_t i = 0, n = geom->getNumGeometries(); i < n; ++i)
      writeLineStringSVG(dynamic_cast<const LineString*>(geom->getGeometryN(i)), out);
  }
  catch (...)
  {
    throw Fmi::Exception::Trace(BCP, "Operation failed!");
  }
}

// ----------------------------------------------------------------------
/*!
 * \brief Handle MultiPolygon
 */
// ----------------------------------------------------------------------

void writeMultiPolygonSVG(const MultiPolygon* geom, std::ostream& out)
{
  try
  {
    if (geom == nullptr || geom->isEmpty())
      return;
    for (size_t i = 0, n = geom->getNumGeometries(); i < n; ++i)
      writePolygonSVG(dynamic_cast<const Polygon*>(geom->getGeometryN(i)), out);
  }
  catch (...)
  {
    throw Fmi::Exception::Trace(BCP, "Operation failed!");
  }
}

// ----------------------------------------------------------------------
/*!
 * \brief Handle GeometryCollection
 */
// ----------------------------------------------------------------------

void writeGeometryCollectionSVG(const GeometryCollection* geom, std::ostream& out)
{
  try
  {
    if (geom == nullptr || geom->isEmpty())
      return;
    for (size_t i = 0, n = geom->getNumGeometries(); i < n; ++i)
      writeSVG(geom->getGeometryN(i), out);
  }
  catch (...)
  {
    throw Fmi::Exception::Trace(BCP, "Operation failed!");
  }
}

// ----------------------------------------------------------------------
/*!
 * \brief Write the geometry as a SVG path
 *
 * This could be more efficiently done if the geometries supported
 * a visitor.
 */
// ----------------------------------------------------------------------

void writeSVG(const Geometry* geom, std::ostream& out)
{
  try
  {
    if (const Point* point = dynamic_cast<const Point*>(geom))
      writePointSVG(point->getCoordinate(), out);
    else if (const LinearRing* lr = dynamic_cast<const LinearRing*>(geom))
      writeLinearRingSVG(lr, out);
    else if (const LineString* ls = dynamic_cast<const LineString*>(geom))
      writeLineStringSVG(ls, out);
    else if (const Polygon* p = dynamic_cast<const Polygon*>(geom))
      writePolygonSVG(p, out);
    else if (const MultiPoint* mp = dynamic_cast<const MultiPoint*>(geom))
      writeMultiPointSVG(mp, out);
    else if (const MultiLineString* ml = dynamic_cast<const MultiLineString*>(geom))
      writeMultiLineStringSVG(ml, out);
    else if (const MultiPolygon* mpg = dynamic_cast<const MultiPolygon*>(geom))
      writeMultiPolygonSVG(mpg, out);
    else if (const GeometryCollection* g = dynamic_cast<const GeometryCollection*>(geom))
      writeGeometryCollectionSVG(g, out);
    else
      throw Fmi::Exception(BCP, "Encountered an unsupported GEOS geometry component");
  }
  catch (...)
  {
    throw Fmi::Exception::Trace(BCP, "Operation failed!");
  }
}
}  // namespace

namespace GeosTools
{
// ----------------------------------------------------------------------
/*!
 * \brief Convert the geometry to a SVG string
 */
// ----------------------------------------------------------------------

std::string getSVG(const Geometry& theGeom, int thePrecision)
{
  try
  {
    // Actual number of decimals to use is automatically selected if
    // no precision is given or it is negative. The code here is modeled
    // after geos::io::WKTWriter

    int decimals = (thePrecision < 0 ? theGeom.getPrecisionModel()->getMaximumSignificantDigits()
                                     : thePrecision);

    std::ostringstream out;
    out << std::setprecision(decimals >= 0 ? decimals : 0);

    writeSVG(&theGeom, out);
    return out.str();
  }
  catch (...)
  {
    throw Fmi::Exception::Trace(BCP, "Operation failed!");
  }
}

}  // namespace GeosTools
}  // namespace SmartMet
