#include "PaddedGrid.h"
#include "NormalGrid.h"
#include "SavitzkyGolay2D.h"

namespace SmartMet
{
namespace Engine
{
namespace Contour
{
PaddedGrid::PaddedGrid(NFmiDataMatrix<float>& theMatrix,
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

  // Expand bbox by one cell of NaN values in all directions for missing value isobands
  const auto box = itsValidCells.bbox();
  itsBBox[0] = box[0] - 1;
  itsBBox[1] = box[1] - 1;
  itsBBox[2] = box[2] + 1;
  itsBBox[3] = box[3] + 1;
}

// Provide wrap-around capability for world data

double PaddedGrid::x(long i, long j) const
{
  if (i <= mini())
    return -itsMaxCoord;
  if (i >= maxi() + 1)
    return itsMaxCoord;

  auto pos = std::min(std::max(0L, j), maxj());
  return itsCoords.x(i, pos);
}

double PaddedGrid::y(long i, long j) const
{
  if (j <= minj())
    return -itsMaxCoord;
  if (j >= maxj() + 1)
    return itsMaxCoord;

  auto pos = std::min(std::max(0L, i), maxi());
  return itsCoords.y(pos, j);
}

float PaddedGrid::operator()(long i, long j) const
{
  if (i <= mini() || i >= maxi() + 1 || j <= minj() || j >= maxj() + 1)
    return std::numeric_limits<double>::quiet_NaN();
  return itsMatrix[i % itsNX][j];
}

bool PaddedGrid::valid(long i, long j) const
{
  if (i <= mini() || i >= maxi() || j <= minj() || j >= maxj())
    return true;  // but full of NaN
  return itsValidCells(i, j);
}

void PaddedGrid::smooth(std::size_t size, std::size_t degree)
{
  NormalGrid normalgrid(itsMatrix, itsCoords, itsValidCells);
  SavitzkyGolay2D::smooth(normalgrid, size, degree);
}

}  // namespace Contour
}  // namespace Engine
}  // namespace SmartMet
