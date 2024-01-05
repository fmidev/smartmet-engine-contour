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

#include <boost/optional.hpp>
#include <gis/BBox.h>
#include <macgyver/DateTime.h>
#include <spine/Parameter.h>
#include <trax/InterpolationType.h>
#include <vector>

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
 public:
  Options() = delete;
  Options(Spine::Parameter theParam,
          const Fmi::DateTime& theTime,
          std::vector<double> theIsoValues);

  Options(Spine::Parameter theParam, const Fmi::DateTime& theTime, std::vector<Range> theLimits);

  std::size_t data_hash_value() const;
  std::size_t filtered_data_hash_value() const;

  friend std::size_t hash_value(const Options& theOptions);

  void transformation(double theMultiplier, double theOffset);
  void transform();

  Trax::InterpolationType interpolation = Trax::InterpolationType::Linear;
  int extrapolation = 0;
  Spine::Parameter parameter;
  Fmi::DateTime time;
  boost::optional<double> level;

  std::vector<double> isovalues;  // for isolines

  std::vector<Range> limits;  // for fills

  boost::optional<Fmi::BBox> bbox;  // user specified WMS file bbox

  boost::optional<double> multiplier;  // for unit conversion
  boost::optional<double> offset;

  boost::optional<std::size_t> filter_size;  // 2D Savitzky-Golay smoother
  boost::optional<std::size_t> filter_degree;

  boost::optional<double> minarea;  // km^2

  bool closed_range = true;  // is last isoband actually 90 <= x <= 100 instead of 90 <= x < 100
  bool validate = false;     // validate the contours according to OGC rules?
  bool strict = false;       // require strict success in contouring, and if validation fails error
  bool desliver = false;     // remove slivers

 private:
  bool hasTransformation() const;

};  // class Options

}  // namespace Contour
}  // namespace Engine
}  // namespace SmartMet
