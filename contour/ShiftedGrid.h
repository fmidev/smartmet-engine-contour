#pragma once

#include <fmt/format.h>
#include <gis/BoolMatrix.h>
#include <gis/CoordinateMatrix.h>
#include <macgyver/Exception.h>
#include <newbase/NFmiDataMatrix.h>
#include <newbase/NFmiPoint.h>
#include <trax/Grid.h>
#include <limits>

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
class ShiftedGrid : public Trax::Grid
{
 public:
  ShiftedGrid() = delete;

  ShiftedGrid(NFmiDataMatrix<float>& theMatrix,
              const Fmi::CoordinateMatrix& theCoords,
              const Fmi::BoolMatrix& theValidCells,
              std::size_t theShift)
      : itsCoords(theCoords),
        itsValidCells(theValidCells),
        itsMatrix(theMatrix),
        itsNX(theMatrix.NX()),
        itsWidth(theCoords.width()),
        itsHeight(theCoords.height()),
        itsShift(theShift)
  {
    if (theCoords.height() != theMatrix.NY() ||
        (theCoords.width() != theMatrix.NX() && theCoords.width() != theMatrix.NX() + 1))
      throw Fmi::Exception(
          BCP,
          fmt::format("Contoured data {}x{} and coordinate dimensions {}x{} mismatch",
                      theMatrix.NX(),
                      theMatrix.NY(),
                      theCoords.width(),
                      theCoords.height()));
    if (theShift == 0)
      throw Fmi::Exception(BCP, "Cannot create shifted grid with shift=0");
  }

  // Provide wrap-around capability for world data

  double x(long i, long j) const override { return itsCoords.x(i, j); }

  double y(long i, long j) const override { return itsCoords.y(i, j); }

  double operator()(long i, long j) const override { return itsMatrix[i % itsNX][j]; }

  // When creating data we do not use shifts
  void set(long i, long j, double z) override { itsMatrix[i][j] = z; }

  bool valid(long i, long j) const override { return itsValidCells(i, j); }

  std::size_t width() const override { return itsWidth; }
  std::size_t height() const override { return itsHeight; }
  std::size_t shift() const override { return itsShift; }

  std::array<long, 4> bbox() const override
  {
    const auto& box = itsValidCells.bbox();
    return {static_cast<long>(box[0]),
            static_cast<long>(box[1]),
            static_cast<long>(box[2]),
            static_cast<long>(box[3])};
  }

 private:
  const Fmi::CoordinateMatrix& itsCoords;
  const Fmi::BoolMatrix& itsValidCells;
  NFmiDataMatrix<float>& itsMatrix;
  const std::size_t itsNX;     // coordinates width
  const std::size_t itsWidth;  // data width
  const std::size_t itsHeight;
  const std::size_t itsShift;  // horizontal shift in start position

};  // class Grid
}  // namespace Contour
}  // namespace Engine
}  // namespace SmartMet
