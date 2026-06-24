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
  ~Config() = default;
  explicit Config(const std::string& theFilename);

  Config() = delete;
  Config(const Config& other) = delete;
  Config& operator=(const Config& other) = delete;
  Config(Config&& other) = delete;
  Config& operator=(Config&& other) = delete;

  int getMaxContourCacheSize() const { return itsMaxContourCacheSize; }

  // Default number of row-bands for parallel contouring (0/1 = single-threaded). The effective
  // value is capped to the number of cores by the engine. Configured via "contour.threads".
  int getThreads() const { return itsThreads; }

 private:
  libconfig::Config itsConfig;
  int itsMaxContourCacheSize;
  int itsThreads = 0;
};

}  // namespace Contour
}  // namespace Engine
}  // namespace SmartMet
