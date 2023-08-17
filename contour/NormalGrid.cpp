#include "NormalGrid.h"
#include "SavitzkyGolay2D.h"

namespace SmartMet
{
namespace Engine
{
namespace Contour
{
NormalGrid::NormalGrid(NFmiDataMatrix<float>& theMatrix,
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
  {
    throw Fmi::Exception(
        BCP,
        fmt::format("Contoured data {}x{} and coordinate dimensions {}x{} mismatch",
                    theMatrix.NX(),
                    theMatrix.NY(),
                    theCoords.width(),
                    theCoords.height()));
  }

  // No expansion as in PaddedGrid
  const auto box = itsValidCells.bbox();
  itsBBox[0] = box[0];
  itsBBox[1] = box[1];
  itsBBox[2] = box[2];
  itsBBox[3] = box[3];
}

// Provide wrap-around capability for world data

double NormalGrid::x(long i, long j) const
{
  if (i < 0 || static_cast<std::size_t>(i) >= itsCoords.width())
    throw Fmi::Exception(
        BCP, fmt::format("NormalGrid accessing coordinates {},{} outside original grid", i, j));

  return itsCoords.x(i, j);
}

double NormalGrid::y(long i, long j) const
{
  if (j < 0 || static_cast<std::size_t>(j) >= itsCoords.height())
    throw Fmi::Exception(
        BCP, fmt::format("NormalGrid accessing coordinates {},{} outside original grid", i, j));

  return itsCoords.y(i, j);
}

float NormalGrid::operator()(long i, long j) const
{
  // Note > instead of >= for max i to provide global wraparound capability
  if (i < 0 || static_cast<std::size_t>(i) >= itsMatrix.NX() || j < 0 ||
      static_cast<std::size_t>(j) >= itsMatrix.NY())
    throw Fmi::Exception(
        BCP, fmt::format("NormalGrid accessing coordinates {},{} outside original grid", i, j));

  return itsMatrix[i % itsNX][j];
}

void NormalGrid::smooth(std::size_t size, std::size_t degree)
{
  NormalGrid normalgrid(itsMatrix, itsCoords, itsValidCells);
  SavitzkyGolay2D::smooth(normalgrid, size, degree);
}

}  // namespace Contour
}  // namespace Engine
}  // namespace SmartMet
