#define BOOST_TEST_MODULE InstantaneousFrequencyObserver
#include <boost/test/included/unit_test.hpp>
#include <cmath>
#include <memory>
#include <vector>
#include <nonlinear_systems/systems/generic_system.hpp>
#include <nonlinear_systems/odes/harmonic_oscillator_ode.hpp>
#include <nonlinear_systems/observers/instantaneous_frequency_observer.hpp>

typedef std::vector<double> state_type;

using namespace nonlinear_systems;

struct Harmonic {
  Harmonic() {
    omega = 2.;
    N = 1;
    d = 2;
    system = std::unique_ptr<GenericSystem<HarmonicOscillatorODE> >(
        new GenericSystem<HarmonicOscillatorODE>(N, d, &omega));
    state_type x_0 = {1., 0.};
    system->SetPosition(x_0);
    dt = 0.01;
  }

  double omega, dt;
  unsigned int N, d;
  std::unique_ptr<GenericSystem<HarmonicOscillatorODE> > system;
  state_type x_0;
  std::vector<state_type> instant_frequency;
  std::vector<double> t;
};


// TODO: Test for multiple oscillators (using uncoupled Kuramoto?)
BOOST_FIXTURE_TEST_CASE(instantaneous_frequency_observer, Harmonic) {
  // for an harmonic oscillator with the initial condition (1, 0) we get the
  // solution x(t) = cos(omega*t), x'(t) = -omega*sin(omega*t)
  // And a solution for the instantaneous frequency f
  // f = -omega^2/(1+(omega^2-1)*sin(omega*t)^2)
  // we start at time dt, because we need the previous point t=0 for the
  // calculation of the derivative using the two-point method
  state_type analytical = {-3.9952, -3.9809, -3.9573, -3.9248};
  unsigned int n_step = 3;
  system->Integrate(dt, n_step, InstantaneousFrequencyObserver<
      GenericSystem<HarmonicOscillatorODE>, state_type>((*system), dt, d, 
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


BOOST_FIXTURE_TEST_CASE(instantaneous_frequency_observer_manual_constructor, 
    Harmonic) {
  // for an harmonic oscillator with the initial condition (1, 0) we get the
  // solution x(t) = cos(omega*t), x'(t) = -omega*sin(omega*t)
  // And a solution for the instantaneous frequency f
  // f = -omega^2/(1+(omega^2-1)*sin(omega*t)^2)
  // we start at time dt, because we need the previous point t=0 for the
  // calculation of the derivative using the two-point method
  state_type analytical = {-3.9952, -3.9809, -3.9573, -3.9248};
  unsigned int n_step = 3;
  state_type system_two_before = system->GetPositionSpherical();
  system_two_before.erase(system_two_before.begin());
  system->Integrate(dt, 1);
  state_type system_one_before = system->GetPositionSpherical();
  system_one_before.erase(system_one_before.begin());
  double t_before = system->GetTime();
  system->Integrate(dt, 1);
  system->Integrate(dt, n_step, InstantaneousFrequencyObserver<
      GenericSystem<HarmonicOscillatorODE>, state_type>((*system), 
        system_one_before, t_before, system_two_before, dt, d, 
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
