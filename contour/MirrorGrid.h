// ======================================================================
/*!
 * \brief Matrix adapter with mirror boundary conditions
 *
 * The data is mirrored at the borders so that the trend in the
 * data is preserved.
 *
 * For example, in 1D case we interpret like this:
 *
 *   f(-2) = f(0) - (f(2)-f(0)) = 2*f(0)-f(2)
 *
 * Similarly at the other end we have
 *
 *   f(w) = f(w-1) - (f(w-2) - f(w-1)) = 2*f(w-1) - f(w-2)
 *
 * To generalize:
 *
 *   f(j) = 2*f(w-1) - f( (w-1)-(j-(w-1))) = 2*f(w-1) - f(2*w-j-2)
 *
 * In the 2D case we simply apply the formulas first for i, then
 * for j.
 *
 * Warning: The mirroring does not extend beyond one grid width!
 */
// ======================================================================

#pragma once

#include "NormalGrid.h"

namespace SmartMet
{
namespace Engine
{
namespace Contour
{
class MirrorGrid
{
 public:
  explicit MirrorGrid(const NormalGrid& theMatrix)
      : M(theMatrix), W(theMatrix.width()), H(theMatrix.height())
  {
  }
  std::size_t width() const { return W; }
  std::size_t height() const { return H; }
  float operator()(long i, long j) const;

 private:
  const NormalGrid& M;
  long W;
  long H;
};

}  // namespace Contour
}  // namespace Engine
}  // namespace SmartMet
