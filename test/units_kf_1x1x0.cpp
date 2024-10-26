#include "fcarouge/kalman.hpp"
#include "fcarouge/unit.hpp"

#include <cassert>
#include <cmath>

namespace fcarouge::test {
namespace {
//! @test Verifies compatibility with the `mp-units` quantities and units
//! library for C++
[[maybe_unused]] auto sample{[] {
  kalman filter{state{60. * m}, output<decltype(0. * m)>,
                estimate_uncertainty{225. * m * m},
                output_uncertainty{25. * m * m}};

  assert(60 * m == filter.x() &&
         "Since our system's dynamic model is constant, i.e. the building "
         "doesn't change its height: 60 meters.");
  assert(225 * m * m == filter.p() &&
         "The extrapolated estimate uncertainty (variance) also doesn't "
         "change: 225");

  filter.update(48.54 * m);

  filter.update(47.11 * m);
  filter.update(55.01 * m);
  filter.update(55.15 * m);
  filter.update(49.89 * m);
  filter.update(40.85 * m);
  filter.update(46.72 * m);
  filter.update(50.05 * m);
  filter.update(51.27 * m);
  filter.update(49.95 * m);

  assert(abs(1 - filter.x() / (49.57 * m)) < 0.001 &&
         "After 10 measurement and update iterations, the building estimated "
         "height is: 49.57m.");

  return 0;
}()};
} // namespace
} // namespace fcarouge::test
