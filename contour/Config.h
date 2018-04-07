#pragma once

#include <boost/utility.hpp>
#include <libconfig.h++>

namespace SmartMet
{
namespace Engine
{
namespace Contour
{
class Config : private boost::noncopyable
{
 public:
  Config() = delete;
  Config(const std::string& theFilename);

  int getMaxContourCacheSize() const { return itsMaxContourCacheSize; }

 private:
  libconfig::Config itsConfig;
  int itsMaxContourCacheSize;
};

}  // namespace Contour
}  // namespace Engine
}  // namespace SmartMet
