#pragma once

#include <libconfig.h++>

namespace SmartMet
{
namespace Engine
{
namespace Contour
{
class Config
{
 public:
  Config() = delete;
  Config(const Config& other) = delete;
  Config& operator=(const Config& other) = delete;
  Config(const std::string& theFilename);

  int getMaxContourCacheSize() const { return itsMaxContourCacheSize; }

 private:
  libconfig::Config itsConfig;
  int itsMaxContourCacheSize;
};

}  // namespace Contour
}  // namespace Engine
}  // namespace SmartMet
