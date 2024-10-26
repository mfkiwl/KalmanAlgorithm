#include "fcarouge/kalman.hpp"
#include "fcarouge/physical_linalg.hpp"

#include <cassert>
#include <cmath>

///////////////////////////////////////////////////////////////////////////////
// https://meetingcpp.com/mcpp/slides/2021/Physical-units-for-matrices6397.pdf
// https://www.youtube.com/watch?v=5dhFtSu3wCo

namespace fcarouge {
template <typename Quantity> struct index {
  Quantity q;

  [[nodiscard]] inline constexpr operator double() const { return 0.; }
};

template <typename... RowIndexes>
using vector =
    physical_column_vector<column_vector<double, sizeof...(RowIndexes)>,
                           RowIndexes...>;

namespace {
[[maybe_unused]] auto test{[] {
  using distance = index<decltype(0. * m)>;
  using velocity = index<decltype(0. * m / s)>;
  using acceleration = index<decltype(0. * m / s2)>;

  using state = vector<distance, velocity, acceleration, distance, velocity,
                       acceleration>;

  using output = vector<distance, distance>;
  
  using estimate_uncertainty = deduce_matrix<state, state>;
  using process_uncertainty = estimate_uncertainty;
  using output_uncertainty = deduce_matrix<output, output>;
  using output_model = deduce_matrix<output, state>;
  using state_transition = state;

  index px0{0. * m};
  index vx0{0. * m / s};
  index ax0{0. * m / s2};
  index py0{0. * m};
  index vy0{0. * m / s};
  index ay0{0. * m / s2};

  state x0{px0, vx0, ax0, py0, vy0, ay0};

  // estimate_uncertainty p0{estimate_uncertainty::zero()};
  // process_uncertainty q0{identity_v<process_uncertainty>};
  // output_uncertainty r0{identity_v<output_uncertainty>};
  // output_model h0{identity_v<output_model>};
  // state_transition f0{identity_v<state_transition>};

  // kalman filter{state{x0},
  //               output<output>,
  //               estimate_uncertainty{p0},
  //               process_uncertainty{q0},
  //               output_uncertainty{r0},
  //               output_model{h0},
  //               state_transition{f0}};

  return 0;
}()};
} // namespace
} // namespace fcarouge

///////////////////////////////////////////////////////////////////////////////
/*
template <typename Numerator, fcarouge::algebraic Denominator>
constexpr auto fcarouge::operator/(const Numerator &lhs, const Denominator &rhs)
    -> fcarouge::deduce_matrix<Numerator, Denominator> {
  return rhs.transpose()
      .fullPivHouseholderQr()
      .solve(lhs.transpose())
      .transpose();
}

namespace fcarouge::sample {
namespace {
template <auto Size> using vector = column_vector<double, Size>;
template <auto Row, auto Column> using matrix = matrix<double, Row, Column>;
using state = fcarouge::state<vector<6>>;

[[maybe_unused]] auto sample{[] {
  using mp_units::si::unit_symbols::m;
  using mp_units::si::unit_symbols::s;

  // A 6x2x0 filter, constant acceleration dynamic model, no control.
  kalman filter{
      // The state X is chosen to be the position, velocity, acceleration in the
      // XY plane: [px, vx, ax, py, vy, ay]. We don't know the vehicle location;
      // we will set initial position, velocity and acceleration to 0.
      state{0., 0., 0., 0., 0., 0.},
      // The vehicle has an onboard location sensor that reports output Z as X
      // and Y
      // coordinates of the system.
      output<vector<2>>,
      // The estimate uncertainty matrix P.
      // Since our initial state vector is a guess, we will set a very high
      // estimate uncertainty. The high estimate uncertainty results in a high
      // Kalman Gain, giving a high weight to the measurement.
      estimate_uncertainty{{500., 0., 0., 0., 0., 0.},
                           {0., 500., 0., 0., 0., 0.},
                           {0., 0., 500., 0., 0., 0.},
                           {0., 0., 0., 500., 0., 0.},
                           {0., 0., 0., 0., 500., 0.},
                           {0., 0., 0., 0., 0., 500.}},
      // The process uncertainty noise matrix Q, constant, computed in place,
      // with  random acceleration standard deviation: σa = 0.2 m.s^-2.
      process_uncertainty{0.2 * 0.2 *
                          matrix<6, 6>{{0.25, 0.5, 0.5, 0., 0., 0.},
                                       {0.5, 1., 1., 0., 0., 0.},
                                       {0.5, 1., 1., 0., 0., 0.},
                                       {0., 0., 0., 0.25, 0.5, 0.5},
                                       {0., 0., 0., 0.5, 1., 1.},
                                       {0., 0., 0., 0.5, 1., 1.}}},
      // The output uncertainty matrix R. Assume that the x and y measurements
      // are uncorrelated, i.e. error in the x coordinate measurement doesn't
      // depend on the error in the y coordinate measurement. In real-life
      // applications, the measurement uncertainty can differ between
      // measurements. In many systems the measurement uncertainty depends on
      // the measurement SNR (signal-to-noise ratio), angle between sensor (or
      // sensors) and target, signal frequency and many other parameters.
      // For the sake of the example simplicity, we will assume a constant
      // measurement uncertainty: R1 = R2...Rn-1 = Rn = R The measurement error
      // standard deviation: σxm = σym = 3m. The variance 9.
      output_uncertainty{{9., 0.}, {0., 9.}},
      // The output model matrix H. The dimension of zn is 2x1 and the dimension
      // of xn is 6x1. Therefore the dimension of the observation matrix H shall
      // be 2x6.
      output_model{{1., 0., 0., 0., 0., 0.}, {0., 0., 0., 1., 0., 0.}},
      // // The state transition matrix F would be:
      state_transition{{1., 1., 0.5, 0., 0., 0.},
                       {0., 1., 1., 0., 0., 0.},
                       {0., 0., 1., 0., 0., 0.},
                       {0., 0., 0., 1., 1., 0.5},
                       {0., 0., 0., 0., 1., 1.},
                       {0., 0., 0., 0., 0., 1.}}};

  // Now we can predict the next state based on the initialization values.
  filter.predict();

  // The measurement values: z1 = [-393.66, 300.4]
  filter.update(-393.66, 300.4);
  filter.predict();

  // And so on, run a step of the filter, predicting and updating, every
  // measurements period: Δt = 1s (constant, built-in).
  const auto step{[&filter](double position_x, double position_y) {
    filter.update(position_x, position_y);
    filter.predict();
  }};

  step(-375.93, 301.78);

  // Verify the example estimated state at 0.1% accuracy.
  assert(std::abs(1 - filter.x()[0] / -277.8) < 0.001 &&
         std::abs(1 - filter.x()[1] / 148.3) < 0.001 &&
         std::abs(1 - filter.x()[2] / 94.5) < 0.001 &&
         std::abs(1 - filter.x()[3] / 249.8) < 0.001 &&
         std::abs(1 - filter.x()[4] / -85.9) < 0.001 &&
         std::abs(1 - filter.x()[5] / -63.6) < 0.001 &&
         "The state estimates expected at 0.1% accuracy.");

  step(-351.04, 295.1);
  step(-328.96, 305.19);
  step(-299.35, 301.06);
  step(-273.36, 302.05);
  step(-245.89, 300);
  step(-222.58, 303.57);
  step(-198.03, 296.33);
  step(-174.17, 297.65);
  step(-146.32, 297.41);
  step(-123.72, 299.61);
  step(-103.47, 299.6);
  step(-78.23, 302.39);
  step(-52.63, 295.04);
  step(-23.34, 300.09);
  step(25.96, 294.72);
  step(49.72, 298.61);
  step(76.94, 294.64);
  step(95.38, 284.88);
  step(119.83, 272.82);
  step(144.01, 264.93);
  step(161.84, 251.46);
  step(180.56, 241.27);
  step(201.42, 222.98);
  step(222.62, 203.73);
  step(239.4, 184.1);
  step(252.51, 166.12);
  step(266.26, 138.71);
  step(271.75, 119.71);
  step(277.4, 100.41);
  step(294.12, 79.76);
  step(301.23, 50.62);
  step(291.8, 32.99);
  step(299.89, 2.14);

  assert(std::abs(1 - filter.x()[0] / 298.5) < 0.006 &&
         std::abs(1 - filter.x()[1] / -1.65) < 0.006 &&
         std::abs(1 - filter.x()[2] / -1.9) < 0.006 &&
         std::abs(1 - filter.x()[3] / -22.5) < 0.006 &&
         std::abs(1 - filter.x()[4] / -26.1) < 0.006 &&
         std::abs(1 - filter.x()[5] / -0.64) < 0.006 &&
         "The state estimates expected at 0.6% accuracy.");
  assert(std::abs(1 - filter.p()(0, 0) / 11.25) < 0.001 &&
         std::abs(1 - filter.p()(0, 1) / 4.5) < 0.001 &&
         std::abs(1 - filter.p()(0, 2) / 0.9) < 0.001 &&
         std::abs(1 - filter.p()(1, 1) / 2.4) < 0.001 &&
         std::abs(1 - filter.p()(2, 2) / 0.2) < 0.001 &&
         std::abs(1 - filter.p()(3, 3) / 11.25) < 0.001 &&
         "The estimate uncertainty expected at 0.1% accuracy."
         "At this point, the position uncertainty px = py = 5, which means "
         "that the standard deviation of the prediction is square root of 5m.");

  return 0;
}()};
} // namespace
} // namespace fcarouge::sample
*/
