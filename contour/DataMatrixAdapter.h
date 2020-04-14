#pragma once

#include <gis/BoolMatrix.h>
#include <gis/CoordinateMatrix.h>
#include <newbase/NFmiDataMatrix.h>
#include <newbase/NFmiPoint.h>

class DataMatrixAdapter
{
 public:
  using value_type = float;
  using coord_type = double;
  using size_type = NFmiDataMatrix<float>::size_type;

  DataMatrixAdapter(NFmiDataMatrix<float>& theMatrix,
                    const Fmi::CoordinateMatrix& theCoords,
                    const Fmi::BoolMatrix& theValidCells)
      : itsCoords(theCoords),
        itsValidCells(theValidCells),
        itsMatrix(theMatrix),
        itsWidth(theMatrix.NX()),
        itsHeight(theMatrix.NY())
  {
  }

  // Provide wrap-around capability for world data
  const value_type& operator()(size_type i, size_type j) const
  {
    return itsMatrix[i % itsWidth][j];
  }

  // Provide wrap-around capability for world data
  value_type& operator()(size_type i, size_type j) { return itsMatrix[i % itsWidth][j]; }
  // Now wrap-around for coordinates, we need both left and right
  // edge coordinates for world data
  coord_type x(size_type i, size_type j) const
  {
    if (i < itsWidth)
      return itsCoords.x(i, j);
    else
      return 360;  // TODO: Could be 180 too for some data
  }

  // For latitude wrap-around value should be OK or the data is not OK
  coord_type y(size_type i, size_type j) const { return itsCoords.y(i % itsWidth, j); }

  bool valid(size_type i, size_type j) const { return itsValidCells(i, j); }

  size_type width() const { return itsWidth; }
  size_type height() const { return itsHeight; }

 private:
  DataMatrixAdapter();
  const Fmi::CoordinateMatrix& itsCoords;
  const Fmi::BoolMatrix& itsValidCells;
  NFmiDataMatrix<float>& itsMatrix;
  const size_type itsWidth;
  const size_type itsHeight;

};  // class DataMatrixAdapter
