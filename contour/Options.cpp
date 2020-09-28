// ======================================================================
/*!
 * \brief Contouring options
 */
// ======================================================================

#include "Options.h"
#include <boost/functional/hash.hpp>
#include <macgyver/Exception.h>
#include <spine/Hash.h>

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

Options::Options(const Spine::Parameter& theParam,
                 const boost::posix_time::ptime& theTime,
                 const std::vector<double>& theIsoValues)
    : interpolation(Linear), parameter(theParam), time(theTime), isovalues(theIsoValues)
{
}

// ----------------------------------------------------------------------
/*!
 * \brief Default options for isobands
 */
// ----------------------------------------------------------------------

Options::Options(const Spine::Parameter& theParam,
                 const boost::posix_time::ptime& theTime,
                 const std::vector<Range>& theLimits)
    : interpolation(Linear), parameter(theParam), time(theTime), limits(theLimits)
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
  boost::hash_combine(seed, theOptions.interpolation);
  boost::hash_combine(seed, theOptions.extrapolation);
  boost::hash_combine(seed, theOptions.parameter);
  boost::hash_combine(seed, theOptions.time);
  boost::hash_combine(seed, theOptions.level);

  for (auto isovalue : theOptions.isovalues)
    boost::hash_combine(seed, isovalue);
  for (auto range : theOptions.limits)
  {
    boost::hash_combine(seed, range.lolimit);
    boost::hash_combine(seed, range.hilimit);
  }

  boost::hash_combine(seed, theOptions.multiplier);
  boost::hash_combine(seed, theOptions.offset);
  boost::hash_combine(seed, theOptions.filter_size);
  boost::hash_combine(seed, theOptions.filter_degree);
  boost::hash_combine(seed, theOptions.minarea);
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
    // NO: boost::hash_combine(seed, interpolation);
    // NO: boost::hash_combine(seed, extrapolation);
    boost::hash_combine(seed, parameter);
    boost::hash_combine(seed, time);
    boost::hash_combine(seed, level);
    // NO: boost::hash_combine(seed, isovalue);
    // NO: boost::hash_combine(seed, lolimit);
    // NO: boost::hash_combine(seed, hilimit);
    // NO: boost::hash_combine(seed, multiplier);
    // NO: boost::hash_combine(seed, offset);
    // NO: boost::hash_combine(seed, filter_size);
    // NO: boost::hash_combine(seed, filter_degree);
    // NO: boost::hash_combine(seed, minarea);
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
    boost::hash_combine(seed, interpolation);
    boost::hash_combine(seed, extrapolation);
    boost::hash_combine(seed, parameter);
    boost::hash_combine(seed, time);
    boost::hash_combine(seed, level);
    // NO: boost::hash_combine(seed, isovalue);
    // NO: boost::hash_combine(seed, lolimit);
    // NO: boost::hash_combine(seed, hilimit);
    boost::hash_combine(seed, multiplier);
    boost::hash_combine(seed, offset);
    boost::hash_combine(seed, filter_size);
    boost::hash_combine(seed, filter_degree);
    boost::hash_combine(seed, minarea);
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
