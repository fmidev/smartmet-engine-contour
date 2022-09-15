#pragma once

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
class Grid : public Trax::Grid
{
 public:
  Grid() = delete;

  Grid(NFmiDataMatrix<float>& theMatrix,
       const Fmi::CoordinateMatrix& theCoords,
       const Fmi::BoolMatrix& theValidCells,
       std::size_t theShift = 0)
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
      throw Fmi::Exception(BCP, "Contoured data and coordinate dimensions mismatch");
  }

  void shell(double maxcoord) { itsMaxCoord = maxcoord; }

  long mini() const { return static_cast<long>(itsValidCells.bbox()[0]); }
  long minj() const { return static_cast<long>(itsValidCells.bbox()[1]); }
  long maxi() const { return static_cast<long>(itsValidCells.bbox()[2] + 2); }
  long maxj() const { return static_cast<long>(itsValidCells.bbox()[3] + 2); }

  // Provide wrap-around capability for world data

  double x(long i, long j) const override
  {
    double ret = 666.0;
    if (i <= mini())
      ret = -itsMaxCoord;
    else if (i >= maxi() + 1)
      ret = itsMaxCoord;
    else
      ret = itsCoords.x(i - 1, std::min(std::max(1L, j), maxj() - 1) - 1);
    return ret;
  }

  double y(long i, long j) const override
  {
    double ret = 666.0;
    if (j <= minj())
      ret = -itsMaxCoord;
    else if (j >= maxj() + 1)
      ret = itsMaxCoord;
    else
      ret = itsCoords.y(std::min(std::max(1L, i), maxi() - 2), j - 1);
    return ret;
  }

  double operator()(long i, long j) const override
  {
    double ret = 666.0;
    if (i <= mini() || i >= maxi() + 1 || j <= minj() || j >= maxj() + 1)
      ret = std::numeric_limits<double>::quiet_NaN();
    else
      ret = itsMatrix[(i - 1) % itsNX][j - 1];
    return ret;
  }

  // When creating data we do not use shifts
  void set(long i, long j, double z) override { itsMatrix[i][j] = z; }

  bool valid(long i, long j) const override
  {
    bool flag = false;

    if (i <= mini() || i >= maxi() + 1 || j <= minj() || j >= maxj() + 1)
      flag = true;  // but full of NaN
    else
      flag = itsValidCells(i - 1, j - 1);
    return flag;
  }

  std::size_t width() const override { return itsWidth; }
  std::size_t height() const override { return itsHeight; }
  std::size_t shift() const override { return itsShift; }

  // We expand the grid by one cell in all directions to surround the data with missing values
  std::array<std::size_t, 4> bbox() const override
  {
    const auto& box = itsValidCells.bbox();
    return {box[0], box[1], box[2] + 2, box[3] + 2};
  }

 private:
  const Fmi::CoordinateMatrix& itsCoords;
  const Fmi::BoolMatrix& itsValidCells;
  NFmiDataMatrix<float>& itsMatrix;
  const std::size_t itsNX;     // coordinates width
  const std::size_t itsWidth;  // data width
  const std::size_t itsHeight;
  const std::size_t itsShift;  // horizontal shift in start position

  double itsMaxCoord = 1e10;

};  // class Grid
}  // namespace Contour
}  // namespace Engine
}  // namespace SmartMet
