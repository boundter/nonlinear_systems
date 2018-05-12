#define BOOST_TEST_MODULE InstantaneousFrequencyObserver
#include <boost/test/included/unit_test.hpp>
#include <cmath>
#include <vector>
#include <nonlinear_systems/systems/generic_system.hpp>
#include <nonlinear_systems/odes/harmonic_oscillator_ode.hpp>
#include <nonlinear_systems/observers/instantaneous_frequency_observer.hpp>

typedef std::vector<double> state_type;

using namespace nonlinear_systems;

// TODO: Test for multiple oscillators (using uncoupled Kuramoto?)
BOOST_AUTO_TEST_CASE(instantaneous_frequency_observer) {
  // for an harmonic oscillator with the initial condition (1, 0) we get the
  // solution x(t) = cos(omega*t), x'(t) = -omega*sin(omega*t)
  // And a solution for the instantaneous frequency f
  // f = -omega^2/(1+(omega^2-1)*sin(omega*t)^2)
  double omega = 2.;
  unsigned int N = 1;
  unsigned int d = 2;
  GenericSystem<HarmonicOscillatorODE> system(N, d, &omega);
  state_type x_0 = {1., 0.};
  system.SetPosition(x_0);
  double dt = 0.01;
  // we start at time dt, because we need the previous point t=0 for the
  // calculation of the derivative using the two-point method
  state_type analytical = {-3.9952, -3.9809, -3.9573, -3.9248};
  std::vector<state_type> instant_frequency;
  std::vector<double> t;
  unsigned int n_step = 3;
  system.Integrate(dt, n_step, InstantaneousFrequencyObserver<
      GenericSystem<HarmonicOscillatorODE>, state_type>(system, dt, d, 
        instant_frequency, t));
  BOOST_REQUIRE_EQUAL(instant_frequency.size(), n_step+1);
  BOOST_REQUIRE_EQUAL(t.size(), n_step+1);
  BOOST_REQUIRE_EQUAL(instant_frequency[0].size(), 1);
  BOOST_CHECK_CLOSE(instant_frequency[0][0], analytical[0], 0.1);
  BOOST_CHECK_CLOSE(t[0], dt, 0.01);
  BOOST_REQUIRE_EQUAL(instant_frequency[1].size(), 1);
  BOOST_CHECK_CLOSE(instant_frequency[1][0], analytical[1], 0.1);
  BOOST_CHECK_CLOSE(t[1], 2*dt, 0.01);
  BOOST_REQUIRE_EQUAL(instant_frequency[2].size(), 1);
  BOOST_CHECK_CLOSE(instant_frequency[2][0], analytical[2], 0.1);
  BOOST_CHECK_CLOSE(t[2], 3*dt, 0.01);
  BOOST_REQUIRE_EQUAL(instant_frequency[3].size(), 1);
  BOOST_CHECK_CLOSE(instant_frequency[3][0], analytical[3], 0.1);
  BOOST_CHECK_CLOSE(t[3], 4*dt, 0.01);
}
