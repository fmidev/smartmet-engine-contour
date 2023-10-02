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
  ~BaseGrid() override;
  BaseGrid() = default;

  BaseGrid(const BaseGrid &) = delete;
  BaseGrid(BaseGrid &&) = delete;
  BaseGrid &operator=(const BaseGrid &) = delete;
  BaseGrid &operator=(BaseGrid &&) = delete;

  virtual void smooth(std::size_t size, std::size_t degree) = 0;
};

}  // namespace Contour
}  // namespace Engine
}  // namespace SmartMet
