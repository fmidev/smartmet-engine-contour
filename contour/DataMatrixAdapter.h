#pragma once

#include <gis/BoolMatrix.h>
#include <gis/CoordinateMatrix.h>
#include <newbase/NFmiDataMatrix.h>
#include <newbase/NFmiPoint.h>
#include <macgyver/Exception.h>

// Note: The values may be one column narrower than the coordinates if the data needs a wraparound
// to fill the globe. In this case the last coordinate is the projected metric coordinate with
// different value to column zero, but essentially the same location on earth in geographic
// coordinates. The values have not been wrapped similarly, hence an extra test is needed when
// fetching values.

namespace SmartMet
{
namespace Engine
{
namespace Contour
{
class DataMatrixAdapter
{
 public:
  using value_type = float;
  using coord_type = double;
  using size_type = NFmiDataMatrix<float>::size_type;

  DataMatrixAdapter(NFmiDataMatrix<float>& theMatrix,
                    const Fmi::CoordinateMatrix& theCoords,
                    const Fmi::BoolMatrix& theValidCells)
      : itsCoords(theCoords),
        itsValidCells(theValidCells),
        itsMatrix(theMatrix),
        itsNX(theMatrix.NX()),
        itsWidth(theCoords.width()),
        itsHeight(theCoords.height())
  {
    if (theCoords.height() != theMatrix.NY() ||
        (theCoords.width() != theMatrix.NX() && theCoords.width() != theMatrix.NX() + 1))
      throw Fmi::Exception(BCP, "Contoured data and coordinate dimensions mismatch");
  }

  // Provide wrap-around capability for world data

  const value_type& operator()(size_type i, size_type j) const { return itsMatrix[i % itsNX][j]; }
  value_type& operator()(size_type i, size_type j) { return itsMatrix[i % itsNX][j]; }

  coord_type x(size_type i, size_type j) const { return itsCoords.x(i, j); }
  coord_type y(size_type i, size_type j) const { return itsCoords.y(i, j); }

  bool valid(size_type i, size_type j) const { return itsValidCells(i, j); }

  size_type width() const { return itsWidth; }
  size_type height() const { return itsHeight; }

 private:
  DataMatrixAdapter();
  const Fmi::CoordinateMatrix& itsCoords;
  const Fmi::BoolMatrix& itsValidCells;
  NFmiDataMatrix<float>& itsMatrix;
  const size_type itsNX;
  const size_type itsWidth;
  const size_type itsHeight;

};  // class DataMatrixAdapter
}  // namespace Contour
}  // namespace Engine
}  // namespace SmartMet
