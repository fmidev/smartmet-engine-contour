#include "SavitzkyGolay2D.h"
#include "MirrorGrid.h"
#include "SavitzkyGolay2DCoefficients.h"
#include <cmath>
#include <stdexcept>
#include <vector>

namespace SmartMet
{
namespace Engine
{
namespace Contour
{
namespace SavitzkyGolay2D
{
void smooth(NormalGrid& input, std::size_t size, std::size_t degree)
{
  if (size == 0 || degree == 0)
    return;

  if (size > 6)
    size = 6;
  if (degree > 5)
    degree = 5;

  int* factor = SavitzkyGolay2DCoefficients::coeffs[size - 1][degree - 1];
  if (factor == 0)
    return;

  // Reflect data at the borders for better results
  MirrorGrid mirror(input);

  int n = 2 * size + 1;

  int denom = SavitzkyGolay2DCoefficients::denoms[size - 1][degree - 1];

  // Process only the desired part, same algorithm as in Impl.cpp for contouring

  auto bbox = input.bbox();
  auto imin = bbox[0];
  auto jmin = bbox[1];
  auto imax = bbox[2] + 1;
  auto jmax = bbox[3] + 1;

  if (imax < imin || jmax < jmin)
    return;

  auto nx = imax - imin + 1;
  auto ny = jmax - jmin + 1;

  // Determine whether we need to split rows into two parts due to shifted global data.

  long shift = input.shift();
  const auto imid = shift - imin;
  const auto needs_two_passes = (shift != 0 && shift > imin && shift < imax);

  std::vector<long> loop_limits;
  if (!needs_two_passes)
    loop_limits = {imin, imax};
  else
    loop_limits = {imid, nx, 0, imid - 1};

  // Holder for temporary results
  std::vector<float> values;
  values.reserve(nx * ny);

  // First calculate smoothened values
  for (auto jj = jmin; jj <= jmax; jj++)
    for (auto k = 0UL; k < loop_limits.size(); k += 2)
      for (auto ii = loop_limits[k]; ii <= loop_limits[k + 1]; ii++)
      {
        if (std::isnan(input(ii, jj)))
          values.push_back(input(ii, jj));
        else
        {
          float sum = 0;
          int pos = 0;
          for (int j = 0; j < n; j++)
            for (int i = 0; i < n; i++)
              sum += (factor[pos++] * mirror(ii + i - size, jj + j - size));
          if (std::isnan(sum))
            values.push_back(input(ii, jj));  // accept only valid smoothed values
          else
            values.push_back(sum / denom);
        }
      }

  // Then replace old values with new ones
  auto in_pos = 0UL;
  for (auto jj = jmin; jj <= jmax; jj++)
    for (auto k = 0UL; k < loop_limits.size(); k += 2)
      for (auto ii = loop_limits[k]; ii <= loop_limits[k + 1]; ii++)
        input.set(ii, jj, values[in_pos++]);
}
}  // namespace SavitzkyGolay2D
}  // namespace Contour
}  // namespace Engine
}  // namespace SmartMet
