#pragma once

#include <cstddef>
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

  // Capacity of the cache of contoured cell masks in bytes. The masks are one bit per grid cell,
  // and the size of the cache is hence configured in megabytes via "cache.max_valid_cells_mbytes".
  // Note that the capacity is divided evenly between the shards of the cache, so a mask larger
  // than 1/16 of the configured size will not be cached at all.
  std::size_t getMaxValidCellsCacheSize() const
  {
    return static_cast<std::size_t>(itsMaxValidCellsCacheMBytes) * 1024 * 1024;
  }

  // Default number of row-bands for parallel contouring (0/1 = single-threaded). The effective
  // value is capped to the number of cores by the engine. Configured via "contour.threads".
  int getThreads() const { return itsThreads; }

 private:
  libconfig::Config itsConfig;
  int itsMaxContourCacheSize;
  int itsMaxValidCellsCacheMBytes = 256;
  int itsThreads = 0;
};

}  // namespace Contour
}  // namespace Engine
}  // namespace SmartMet
