#pragma once

#include <newbase/NFmiDataMatrix.h>

class DataMatrixAdapter
{
 public:
  typedef float value_type;
  typedef double coord_type;

  typedef NFmiDataMatrix<float>::size_type size_type;

  DataMatrixAdapter(NFmiDataMatrix<float>& theMatrix, const NFmiDataMatrix<NFmiPoint>& theCoords)
      : itsCoords(theCoords),
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
      return itsCoords[i][j].X();
    else
      return 360;  // TODO: Could be 180 too for some data
  }

  // For latitude wrap-around value should be OK or the data is not OK
  coord_type y(size_type i, size_type j) const { return itsCoords[i % itsWidth][j].Y(); }
  size_type width() const { return itsWidth; }
  size_type height() const { return itsHeight; }
  void swap(DataMatrixAdapter& theOther)
  {
    std::swap(itsMatrix, theOther.itsMatrix);
    std::swap(itsWidth, theOther.itsWidth);
    std::swap(itsHeight, theOther.itsHeight);
  }

 private:
  DataMatrixAdapter();
  const NFmiDataMatrix<NFmiPoint>& itsCoords;
  NFmiDataMatrix<float>& itsMatrix;
  size_type itsWidth;
  size_type itsHeight;

};  // class DataMatrixAdapter
