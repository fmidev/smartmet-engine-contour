#pragma once

#include "BaseGrid.h"
#include <fmt/format.h>
#include <gis/BoolMatrix.h>
#include <gis/CoordinateMatrix.h>
#include <macgyver/Exception.h>
#include <newbase/NFmiDataMatrix.h>
#include <newbase/NFmiPoint.h>
#include <limits>

namespace SmartMet
{
namespace Engine
{
namespace Contour
{
// A Grid which has been surrounded by very large cells containing NaN
// values outside the original grid. This enables us to get missing
// data isobands from areas outside the original grid.

class PaddedGrid : public BaseGrid
{
 public:
  virtual ~PaddedGrid() = default;
  PaddedGrid() = delete;

  PaddedGrid(NFmiDataMatrix<float>& theMatrix,
             const Fmi::CoordinateMatrix& theCoords,
             const Fmi::BoolMatrix& theValidCells);

  void shell(double value) { itsMaxCoord = value; }
  double shell() const override { return itsMaxCoord; }

  // Provide wrap-around capability for world data

  double x(long i, long j) const override;
  double y(long i, long j) const override;
  float operator()(long i, long j) const override;
  void set(long i, long j, float z) override { itsMatrix[i][j] = z; }

  bool valid(long i, long j) const override;

  std::size_t width() const override { return itsWidth; }
  std::size_t height() const override { return itsHeight; }

  // We expand the grid by one cell in all directions to surround the data with missing values
  std::array<long, 4> bbox() const override { return itsBBox; }

  void smooth(std::size_t size, std::size_t degree) override;

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

};  // class PaddedGrid
}  // namespace Contour
}  // namespace Engine
}  // namespace SmartMet
