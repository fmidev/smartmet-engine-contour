// ======================================================================
/*!
 * \brief Contouring options
 *
 * 1. If isovalue is set, an isoline will be calculated.
 * 2. Otherwise the range lolimit...hilimit will be calculated.
 *   a) If both lolimit and hilimit are missing, the missing value range
 *      is calculated
 *   b) If only lolimit is missing, the range -infinity...hilimit is used
 *   c) If only hilimit is missing, the range lolimit...infinity is used
 */
// ----------------------------------------------------------------------

#pragma once

#include "Interpolation.h"
#include <spine/Parameter.h>
#include <boost/optional.hpp>
#include <boost/date_time/posix_time/ptime.hpp>

namespace SmartMet
{
namespace Engine
{
namespace Contour
{
struct Range
{
  boost::optional<double> lolimit;
  boost::optional<double> hilimit;

  Range(boost::optional<double> lo, boost::optional<double> hi) : lolimit(lo), hilimit(hi) {}
  Range(double lo, double hi) : lolimit(lo), hilimit(hi) {}
};

struct Options
{
 private:
  Options();

 public:
  Options(const SmartMet::Spine::Parameter& theParam,
          const boost::posix_time::ptime& theTime,
          const std::vector<double>& theIsoValues);

  Options(const SmartMet::Spine::Parameter& theParam,
          const boost::posix_time::ptime& theTime,
          const std::vector<Range>& theLimits);

  std::size_t data_hash_value() const;
  std::size_t filtered_data_hash_value() const;

  friend std::size_t hash_value(const Options& theOptions);

  void transformation(double theMultiplier, double theOffset);
  bool hasTransformation() const;

  Interpolation interpolation;
  SmartMet::Spine::Parameter parameter;
  boost::posix_time::ptime time;
  boost::optional<double> level;

  std::vector<double> isovalues;  // for isolines

  std::vector<Range> limits;  // for fills

  boost::optional<double> multiplier;  // for unit conversion
  boost::optional<double> offset;

  boost::optional<std::size_t> filter_size;  // 2D Savitzky-Golay smoother
  boost::optional<std::size_t> filter_degree;

};  // class Options

}  // namespace Contour
}  // namespace Engine
}  // namespace SmartMet
