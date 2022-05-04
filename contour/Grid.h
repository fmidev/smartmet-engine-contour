#pragma once

#include <gis/BoolMatrix.h>
#include <gis/CoordinateMatrix.h>
#include <macgyver/Exception.h>
#include <newbase/NFmiDataMatrix.h>
#include <newbase/NFmiPoint.h>
#include <trax/Grid.h>

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
class Grid : public Trax::Grid
{
 public:
  Grid() = delete;

  Grid(const NFmiDataMatrix<float>& theMatrix,
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

  double x(long i, long j) const override { return itsCoords.x(i, j); }
  double y(long i, long j) const override { return itsCoords.y(i, j); }
  double operator()(long i, long j) const override { return itsMatrix[i % itsNX][j]; }
  void set(long i, long j, double z) override {}  // not needed

  bool valid(long i, long j) const override
  {
    if (i < 0 || j < 0 || static_cast<std::size_t>(i) >= itsWidth ||
        static_cast<std::size_t>(j) >= itsHeight)
      return false;
    return itsValidCells(i, j);
  }

  std::size_t width() const override { return itsWidth; }
  std::size_t height() const override { return itsHeight; }

  std::array<std::size_t, 4> bbox() const override { return itsValidCells.bbox(); }

 private:
  const Fmi::CoordinateMatrix& itsCoords;
  const Fmi::BoolMatrix& itsValidCells;
  const NFmiDataMatrix<float>& itsMatrix;
  const std::size_t itsNX;     // coordinates width
  const std::size_t itsWidth;  // data width
  const std::size_t itsHeight;

};  // class Grid
}  // namespace Contour
}  // namespace Engine
}  // namespace SmartMet