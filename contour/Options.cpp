// ======================================================================
/*!
 * \brief Contouring options
 */
// ======================================================================

#include "Options.h"
#include <macgyver/Exception.h>
#include <macgyver/Hash.h>

namespace SmartMet
{
namespace Engine
{
namespace Contour
{
// ----------------------------------------------------------------------
/*!
 * \brief Default options for isolines
 */
// ----------------------------------------------------------------------

Options::Options(Spine::Parameter theParam,
                 const boost::posix_time::ptime& theTime,
                 std::vector<double> theIsoValues)
    : parameter(std::move(theParam)), time(theTime), isovalues(std::move(theIsoValues))
{
}

// ----------------------------------------------------------------------
/*!
 * \brief Default options for isobands
 */
// ----------------------------------------------------------------------

Options::Options(Spine::Parameter theParam,
                 const boost::posix_time::ptime& theTime,
                 std::vector<Range> theLimits)
    : parameter(std::move(theParam)), time(theTime), limits(std::move(theLimits))
{
}

// ----------------------------------------------------------------------
/*!
 * \brief Set a linear transformation (mostly for unit conversions)
 */
// ----------------------------------------------------------------------

void Options::transformation(double theMultiplier, double theOffset)
{
  try
  {
    multiplier = theMultiplier;
    offset = theOffset;
  }
  catch (...)
  {
    throw Fmi::Exception::Trace(BCP, "Operation failed!");
  }
}

// ----------------------------------------------------------------------
/*!
 * \brief Test if a transformation has been set
 */
// ----------------------------------------------------------------------

bool Options::hasTransformation() const
{
  return multiplier || offset;
}
// ----------------------------------------------------------------------
/*!
 * \brief Calculate a hash for the request
 */
// ----------------------------------------------------------------------

std::size_t hash_value(const Options& theOptions)
{
  std::size_t seed = 0;
  Fmi::hash_combine(seed, Fmi::hash_value(static_cast<int>(theOptions.interpolation)));
  Fmi::hash_combine(seed, Fmi::hash_value(theOptions.extrapolation));
  Fmi::hash_combine(seed, theOptions.parameter.hashValue());
  Fmi::hash_combine(seed, Fmi::hash_value(theOptions.time));
  Fmi::hash_combine(seed, Fmi::hash_value(theOptions.level));

  for (auto isovalue : theOptions.isovalues)
    Fmi::hash_combine(seed, Fmi::hash_value(isovalue));
  for (auto range : theOptions.limits)
  {
    Fmi::hash_combine(seed, Fmi::hash_value(range.lolimit));
    Fmi::hash_combine(seed, Fmi::hash_value(range.hilimit));
  }

  Fmi::hash_combine(seed, Fmi::hash_value(theOptions.multiplier));
  Fmi::hash_combine(seed, Fmi::hash_value(theOptions.offset));
  Fmi::hash_combine(seed, Fmi::hash_value(theOptions.filter_size));
  Fmi::hash_combine(seed, Fmi::hash_value(theOptions.filter_degree));
  Fmi::hash_combine(seed, Fmi::hash_value(theOptions.minarea));
  return seed;
}

// ----------------------------------------------------------------------
/*!
 * \brief Calculate a hash for the data, ignoring contour limits
 *
 * This is used for caching information relating to the data being
 * contoured.
 */
// ----------------------------------------------------------------------

std::size_t Options::data_hash_value() const
{
  try
  {
    std::size_t seed = 0;
    // NO: Fmi::hash_combine(seed, Fmi::hash_value(interpolation));
    // NO: Fmi::hash_combine(seed, Fmi::hash_value(extrapolation));
    Fmi::hash_combine(seed, parameter.hashValue());
    Fmi::hash_combine(seed, Fmi::hash_value(time));
    Fmi::hash_combine(seed, Fmi::hash_value(level));
    // NO: Fmi::hash_combine(seed, Fmi::hash_value(isovalue));
    // NO: Fmi::hash_combine(seed, Fmi::hash_value(lolimit));
    // NO: Fmi::hash_combine(seed, Fmi::hash_value(hilimit));
    // NO: Fmi::hash_combine(seed, Fmi::hash_value(multiplier));
    // NO: Fmi::hash_combine(seed, Fmi::hash_value(offset));
    // NO: Fmi::hash_combine(seed, Fmi::hash_value(filter_size));
    // NO: Fmi::hash_combine(seed, Fmi::hash_value(filter_degree));
    // NO: Fmi::hash_combine(seed, Fmi::hash_value(minarea));
    return seed;
  }
  catch (...)
  {
    throw Fmi::Exception::Trace(BCP, "Operation failed!");
  }
}

// ----------------------------------------------------------------------
/*!
 * \brief Calculate a hash for the data, ignoring contour limits
 *
 * This is used for caching information relating to the data being
 * contoured.
 */
// ----------------------------------------------------------------------

std::size_t Options::filtered_data_hash_value() const
{
  try
  {
    std::size_t seed = 0;
    Fmi::hash_combine(seed, Fmi::hash_value(static_cast<int>(interpolation)));
    Fmi::hash_combine(seed, Fmi::hash_value(extrapolation));
    Fmi::hash_combine(seed, parameter.hashValue());
    Fmi::hash_combine(seed, Fmi::hash_value(time));
    Fmi::hash_combine(seed, Fmi::hash_value(level));
    // NO: Fmi::hash_combine(seed, Fmi::hash_value(isovalue));
    // NO: Fmi::hash_combine(seed, Fmi::hash_value(lolimit));
    // NO: Fmi::hash_combine(seed, Fmi::hash_value(hilimit));
    Fmi::hash_combine(seed, Fmi::hash_value(multiplier));
    Fmi::hash_combine(seed, Fmi::hash_value(offset));
    Fmi::hash_combine(seed, Fmi::hash_value(filter_size));
    Fmi::hash_combine(seed, Fmi::hash_value(filter_degree));
    Fmi::hash_combine(seed, Fmi::hash_value(minarea));
    return seed;
  }
  catch (...)
  {
    throw Fmi::Exception::Trace(BCP, "Operation failed!");
  }
}

}  // namespace Contour
}  // namespace Engine
}  // namespace SmartMet
