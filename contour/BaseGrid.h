#pragma once

#include <trax/Grid.h>

namespace SmartMet
{
namespace Engine
{
namespace Contour
{
class BaseGrid : public Trax::Grid
{
 public:
  virtual ~BaseGrid() = default;
  BaseGrid() = default;

  virtual void smooth(std::size_t size, std::size_t degree) = 0;
};

}  // namespace Contour
}  // namespace Engine
}  // namespace SmartMet
