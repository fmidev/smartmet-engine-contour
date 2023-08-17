#include "ShiftedGrid.h"
#include "NormalGrid.h"
#include "SavitzkyGolay2D.h"
#include <fmt/format.h>

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
ShiftedGrid::ShiftedGrid(NFmiDataMatrix<float>& theMatrix,
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

  // Convert limits to long as expected by Trax API
  const auto box = itsValidCells.bbox();
  itsBBox[0] = box[0];
  itsBBox[1] = box[1];
  itsBBox[2] = box[2];
  itsBBox[3] = box[3];
}

void ShiftedGrid::smooth(std::size_t size, std::size_t degree)
{
  NormalGrid normalgrid(itsMatrix, itsCoords, itsValidCells);
  SavitzkyGolay2D::smooth(normalgrid, size, degree);
}

}  // namespace Contour
}  // namespace Engine
}  // namespace SmartMet
