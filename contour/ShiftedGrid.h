#pragma once

#include <gis/BoolMatrix.h>
#include <gis/CoordinateMatrix.h>
#include <macgyver/Exception.h>
#include <newbase/NFmiDataMatrix.h>
#include <newbase/NFmiPoint.h>
#include <trax/Grid.h>

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
        itsHeight(theCoords.height()),
        itsShift(theShift)
  {
    // Note that width is taken from data, not coordinates. Shifting usually
    // effectively breaks contouring if we try to provide wraparound capabilities.

    if (theCoords.height() != theMatrix.NY() ||
        (theCoords.width() != theMatrix.NX() && theCoords.width() != theMatrix.NX() + 1))
      throw Fmi::Exception(BCP, "Contoured data and coordinate dimensions mismatch");

    if (theShift == 0)
      throw Fmi::Exception(
          BCP, "A zero shift for ShiftedGrid disables global data wraparound, fix code");
  }

  // Provide wrap-around capability for world data

  double x(long i, long j) const override { return itsCoords.x((i + itsShift) % itsNX, j); }
  double y(long i, long j) const override { return itsCoords.y((i + itsShift) % itsNX, j); }
  double operator()(long i, long j) const override { return itsMatrix[(i + itsShift) % itsNX][j]; }
  void set(long i, long j, double z) override { itsMatrix[i][j] = z; }
  void shift(std::size_t shift) { itsShift = shift; }

  bool valid(long i, long j) const override
  {
    if (i < 0 || j < 0 || static_cast<std::size_t>(i) >= itsNX ||
        static_cast<std::size_t>(j) >= itsHeight)
      return false;
    return itsValidCells((i + itsShift) % itsNX, j);
  }

  std::size_t width() const override { return itsNX; }
  std::size_t height() const override { return itsHeight; }

  std::array<std::size_t, 4> bbox() const override
  {
    auto ret = itsValidCells.bbox();  // xmin,ymin,xmax,ymax
    // Prevent wraparound of cells. One minus for cells vs vertices, another for dropped cell
    ret[2] = std::min(ret[2], itsNX - 2);
    return ret;
  }

 private:
  const Fmi::CoordinateMatrix& itsCoords;
  const Fmi::BoolMatrix& itsValidCells;
  NFmiDataMatrix<float>& itsMatrix;
  const std::size_t itsNX;  // coordinates width
  const std::size_t itsHeight;
  std::size_t itsShift = 0;  // horizontal shift for global data wraparound

};  // class ShiftedGrid
}  // namespace Contour
}  // namespace Engine
}  // namespace SmartMet
