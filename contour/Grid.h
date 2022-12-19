#pragma once

#include <fmt/format.h>
#include <gis/BoolMatrix.h>
#include <gis/CoordinateMatrix.h>
#include <macgyver/Exception.h>
#include <newbase/NFmiDataMatrix.h>
#include <newbase/NFmiPoint.h>
#include <trax/Grid.h>
#include <limits>

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
    const auto bbox = itsValidCells.bbox();
    itsBBox[0] = bbox[0] - 1;
    itsBBox[1] = bbox[1] - 1;
    itsBBox[2] = bbox[2] + 1;
    itsBBox[3] = bbox[3] + 1;
  }

  void shell(double value) { itsMaxCoord = value; }
  double shell() const override { return itsMaxCoord; }

  // Provide wrap-around capability for world data

  double x(long i, long j) const override
  {
    double ret = 666.0;
    if (i <= mini())
      ret = -itsMaxCoord;
    else if (i >= maxi() + 1)
      ret = itsMaxCoord;
    else
    {
      auto pos = std::min(std::max(0L, j), maxj());
      ret = itsCoords.x(i, pos);
    }

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
    {
      auto pos = std::min(std::max(0L, i), maxi());
      ret = itsCoords.y(pos, j);
    }
    return ret;
  }

  double operator()(long i, long j) const override
  {
    double ret = 666.0;
    if (i <= mini() || i >= maxi() + 1 || j <= minj() || j >= maxj() + 1)
      ret = std::numeric_limits<double>::quiet_NaN();
    else
      ret = itsMatrix[i % itsNX][j];
    return ret;
  }

  void set(long i, long j, double z) override { itsMatrix[i][j] = z; }

  bool valid(long i, long j) const override
  {
    bool flag = false;

    if (i <= mini() || i >= maxi() || j <= minj() || j >= maxj())
      flag = true;  // but full of NaN
    else
      flag = itsValidCells(i, j);
    return flag;
  }

  std::size_t width() const override { return itsWidth; }
  std::size_t height() const override { return itsHeight; }

  // We expand the grid by one cell in all directions to surround the data with missing values
  std::array<long, 4> bbox() const override { return itsBBox; }

 private:
  long mini() const { return itsBBox[0]; }
  long minj() const { return itsBBox[1]; }
  long maxi() const { return itsBBox[2]; }
  long maxj() const { return itsBBox[3]; }

  const Fmi::CoordinateMatrix& itsCoords;
  const Fmi::BoolMatrix& itsValidCells;
  NFmiDataMatrix<float>& itsMatrix;
  const std::size_t itsNX;     // coordinates width
  const std::size_t itsWidth;  // data width
  const std::size_t itsHeight;

  std::array<long, 4> itsBBox;  // limits for cell iteration

  double itsMaxCoord = 1e10;

};  // class Grid
}  // namespace Contour
}  // namespace Engine
}  // namespace SmartMet
