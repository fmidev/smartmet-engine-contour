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
  Config(const std::string& theFilename);

  Config() = delete;
  Config(const Config& other) = delete;
  Config& operator=(const Config& other) = delete;
  Config(Config&& other) = delete;
  Config& operator=(Config&& other) = delete;

  int getMaxContourCacheSize() const { return itsMaxContourCacheSize; }

 private:
  libconfig::Config itsConfig;
  int itsMaxContourCacheSize;
};

}  // namespace Contour
}  // namespace Engine
}  // namespace SmartMet
