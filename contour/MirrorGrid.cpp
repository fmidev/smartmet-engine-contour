#include "MirrorGrid.h"
#include <cassert>

namespace SmartMet
{
namespace Engine
{
namespace Contour
{
float MirrorGrid::operator()(long i, long j) const
{
  // Simple reflection does not work out of these bounds
  assert(i > -W);
  assert(i < 2 * W - 1);
  assert(j > -H);
  assert(j < 2 * H - 1);

  if (i < 0)
  {
    if (j < 0)
      return 2 * (2 * M(0, 0) - M(0, -j)) - (2 * M(-i, 0) - M(-i, -j));
    else if (j >= H)
      return 2 * (2 * M(0, H - 1) - M(0, 2 * H - j - 2)) -
             (2 * M(-i, H - 1) - M(-i, 2 * H - j - 2));
    else
      return 2 * M(0, j) - M(-i, j);
  }
  else if (i >= W)
  {
    if (j < 0)
      return 2 * (2 * M(W - 1, 0) - M(W - 1, -j)) -
             (2 * M(2 * W - i - 2, 0) - M(2 * W - i - 2, -j));
    else if (j >= H)
      return 2 * (2 * M(W - 1, H - 1) - M(W - 1, 2 * H - j - 2)) -
             (2 * M(2 * W - i - 2, H - 1) - M(2 * W - i - 2, 2 * H - j - 2));
    else
      return 2 * M(W - 1, j) - M(2 * W - i - 2, j);
  }
  else
  {
    if (j < 0)
      return 2 * M(i, 0) - M(i, -j);
    else if (j >= H)
      return 2 * M(i, H - 1) - M(i, 2 * H - j - 2);
    else
      return M(i, j);
  }
}

}  // namespace Contour
}  // namespace Engine
}  // namespace SmartMet
