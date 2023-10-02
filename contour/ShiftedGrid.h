#pragma once

#include "BaseGrid.h"
#include <gis/BoolMatrix.h>
#include <gis/CoordinateMatrix.h>
#include <macgyver/Exception.h>
#include <newbase/NFmiDataMatrix.h>
#include <newbase/NFmiPoint.h>
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
class ShiftedGrid : public BaseGrid
{
 public:
  ~ShiftedGrid() override = default;
  ShiftedGrid() = delete;

  ShiftedGrid(const ShiftedGrid&) = delete;
  ShiftedGrid(ShiftedGrid&&) = delete;
  ShiftedGrid& operator=(const ShiftedGrid&) = delete;
  ShiftedGrid& operator=(ShiftedGrid&&) = delete;

  ShiftedGrid(NFmiDataMatrix<float>& theMatrix,
              const Fmi::CoordinateMatrix& theCoords,
              const Fmi::BoolMatrix& theValidCells,
              std::size_t theShift);

  // Provide wrap-around capability for world data

  double x(long i, long j) const override { return itsCoords.x(i, j); }

  double y(long i, long j) const override { return itsCoords.y(i, j); }

  float operator()(long i, long j) const override { return itsMatrix[i % itsNX][j]; }

  // When creating data we do not use shifts
  void set(long i, long j, float z) override { itsMatrix[i][j] = z; }

  bool valid(long i, long j) const override { return itsValidCells(i, j); }

  std::size_t width() const override { return itsWidth; }
  std::size_t height() const override { return itsHeight; }
  std::size_t shift() const override { return itsShift; }

  std::array<long, 4> bbox() const override { return itsBBox; }

  void smooth(std::size_t size, std::size_t degree) override;

 private:
  const Fmi::CoordinateMatrix& itsCoords;
  const Fmi::BoolMatrix& itsValidCells;
  NFmiDataMatrix<float>& itsMatrix;
  std::array<long, 4> itsBBox;

  const std::size_t itsNX;     // coordinates width
  const std::size_t itsWidth;  // data width
  const std::size_t itsHeight;
  const std::size_t itsShift;  // horizontal shift in start position

};  // class ShiftedGrid
}  // namespace Contour
}  // namespace Engine
}  // namespace SmartMet
