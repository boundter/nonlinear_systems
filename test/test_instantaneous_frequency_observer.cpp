#define BOOST_TEST_MODULE InstantaneousFrequencyObserver
#include <boost/test/included/unit_test.hpp>
#include <cmath>
#include <memory>
#include <vector>
#include <nonlinear_systems/odes/kuramoto_sakaguchi_ode.hpp>
#include <nonlinear_systems/systems/m_kuramoto_sakaguchi_system.hpp>
#include <nonlinear_systems/systems/generic_system.hpp>
#include <nonlinear_systems/odes/harmonic_oscillator_ode.hpp>
#include <nonlinear_systems/observers/instantaneous_frequency_observer.hpp>

typedef std::vector<double> state_type;

using namespace nonlinear_systems;
// TODO: Refactor in seperate file
class KuramotoSystem: public GenericSystem<KuramotoSakaguchiODE, state_type> {
    public:
      KuramotoSystem(unsigned int N_oscillator, state_type frequency, 
          double coupling, double phase_shift)
      : GenericSystem<KuramotoSakaguchiODE, state_type>(N_oscillator, 1) {
        this->_ode = std::unique_ptr<KuramotoSakaguchiODE>(new 
            KuramotoSakaguchiODE(N_oscillator, frequency, coupling, phase_shift));
      }
};

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


struct Kuramoto {
  Kuramoto() {
    N = 2;
    frequency = {1., 2.};
    analytical_freq = 0.5*(frequency[0] + frequency[1]);
    system = std::unique_ptr<KuramotoSystem>(new KuramotoSystem(N, frequency, 0., 0.));
    system->SetPosition({0., 0.});
    dt = 0.01;
    n = 2;
  }

  unsigned int N;
  state_type frequency;
  double analytical_freq;
  std::unique_ptr<KuramotoSystem> system;
  double dt;
  unsigned int n;
  std::vector<state_type> instant_freq;
  std::vector<state_type> instant_freq_mean_field;
  state_type t;
};


struct MKuramoto {
 MKuramoto() {
    N = {2, 2};
    frequency = {1., 2., 3., 4.};
    analytical_freq = {0.5*(frequency[0] + frequency[1]),
                       0.5*(frequency[2] + frequency[3])};
    coupling = {{0., 0.}, {0., 0.}};
    system = std::unique_ptr<MKuramotoSakaguchiSystem>(
        new MKuramotoSakaguchiSystem(frequency, coupling, coupling, N));
    system->SetPosition({0., 0., 0., 0.});
    dt = 0.01;
    n = 2;
 }
  
 std::vector<unsigned int> N;
 state_type frequency;
 state_type analytical_freq;
 std::vector<state_type> coupling;
 std::unique_ptr<MKuramotoSakaguchiSystem> system;
 double dt;
 unsigned int n;
 std::vector<state_type> instant_freq;
 std::vector<state_type> instant_freq_mean_field;
 state_type t;
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


BOOST_FIXTURE_TEST_CASE(instantaneous_frequency_mean_field_observer, Harmonic) {
  // for an harmonic oscillator with the initial condition (1, 0) we get the
  // solution x(t) = cos(omega*t), x'(t) = -omega*sin(omega*t)
  // And a solution for the instantaneous frequency f
  // f = -omega^2/(1+(omega^2-1)*sin(omega*t)^2)
  // we start at time dt, because we need the previous point t=0 for the
  // calculation of the derivative using the two-point method
  state_type analytical = {-3.9952, -3.9809, -3.9573, -3.9248};
  unsigned int n_step = 3;
  system->Integrate(dt, n_step, InstantaneousFrequencyMeanFieldObserver<
      GenericSystem<HarmonicOscillatorODE>, state_type>((*system), dt, 
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

BOOST_FIXTURE_TEST_CASE(instantaneous_frequency_mean_field_multiple, Kuramoto) {
  system->Integrate(dt, n, InstantaneousFrequencyMeanFieldObserver<KuramotoSystem,
      state_type>((*system), dt, instant_freq, t));
  BOOST_REQUIRE_EQUAL(instant_freq.size(), n+1);
  BOOST_REQUIRE_EQUAL(t.size(), n+1);
  BOOST_REQUIRE_EQUAL(instant_freq[0].size(), 1);
  BOOST_CHECK_CLOSE(instant_freq[0][0], analytical_freq, 0.1);
  BOOST_CHECK_CLOSE(t[0], dt, 0.1);
  BOOST_REQUIRE_EQUAL(instant_freq[1].size(), 1);
  BOOST_CHECK_CLOSE(instant_freq[1][0], analytical_freq, 0.1);
  BOOST_CHECK_CLOSE(t[1], 2*dt, 0.1);
  BOOST_REQUIRE_EQUAL(instant_freq[2].size(), 1);
  BOOST_CHECK_CLOSE(instant_freq[2][0], analytical_freq, 0.1);
  BOOST_CHECK_CLOSE(t[2], 3*dt, 0.1);
}


BOOST_FIXTURE_TEST_CASE(instantaneous_frequency_mean_field_and_phase_multiple,
    Kuramoto) {
  system->Integrate(dt, n, InstantaneousFrequencyMeanFieldAndPhaseObserver<
      KuramotoSystem, state_type>((*system), dt, 1, instant_freq, instant_freq_mean_field,
        t));
  BOOST_REQUIRE_EQUAL(instant_freq.size(), n+1);
  BOOST_REQUIRE_EQUAL(instant_freq_mean_field.size(), n+1);
  BOOST_REQUIRE_EQUAL(t.size(), n+1);
  BOOST_REQUIRE_EQUAL(instant_freq[0].size(), 2);
  BOOST_REQUIRE_EQUAL(instant_freq_mean_field[0].size(), 1);
  BOOST_CHECK_CLOSE(instant_freq[0][0], frequency[0], 0.1);
  BOOST_CHECK_CLOSE(instant_freq[0][1], frequency[1], 0.1);
  BOOST_CHECK_CLOSE(instant_freq_mean_field[0][0], analytical_freq, 0.1);
  BOOST_CHECK_CLOSE(t[0], dt, 0.1);
  BOOST_REQUIRE_EQUAL(instant_freq[1].size(), 2);
  BOOST_REQUIRE_EQUAL(instant_freq_mean_field[1].size(), 1);
  BOOST_CHECK_CLOSE(instant_freq[1][0], frequency[0], 0.1);
  BOOST_CHECK_CLOSE(instant_freq[1][1], frequency[1], 0.1);
  BOOST_CHECK_CLOSE(instant_freq_mean_field[1][0], analytical_freq, 0.1);
  BOOST_CHECK_CLOSE(t[1], 2*dt, 0.1);
  BOOST_REQUIRE_EQUAL(instant_freq[2].size(), 2);
  BOOST_REQUIRE_EQUAL(instant_freq_mean_field[2].size(), 1);
  BOOST_CHECK_CLOSE(instant_freq[2][0], frequency[0], 0.1);
  BOOST_CHECK_CLOSE(instant_freq[2][1], frequency[1], 0.1);
  BOOST_CHECK_CLOSE(instant_freq_mean_field[2][0], analytical_freq, 0.1);
  BOOST_CHECK_CLOSE(t[2], 3*dt, 0.1);
}


BOOST_FIXTURE_TEST_CASE(instantaneous_frequency_mean_field_network, 
    MKuramoto) {
  system->Integrate(dt, n, InstantaneousFrequencyMeanFieldObserver<
      MKuramotoSakaguchiSystem,
      state_type, std::vector<state_type>>((*system), dt, instant_freq, t));
  BOOST_REQUIRE_EQUAL(instant_freq.size(), n+1);
  BOOST_REQUIRE_EQUAL(t.size(), n+1);
  BOOST_REQUIRE_EQUAL(instant_freq[0].size(), 2);
  BOOST_CHECK_CLOSE(instant_freq[0][0], analytical_freq[0], 0.1);
  BOOST_CHECK_CLOSE(instant_freq[0][1], analytical_freq[1], 0.1);
  BOOST_CHECK_CLOSE(t[0], dt, 0.1);
  BOOST_REQUIRE_EQUAL(instant_freq[1].size(), 2);
  BOOST_CHECK_CLOSE(instant_freq[1][0], analytical_freq[0], 0.1);
  BOOST_CHECK_CLOSE(instant_freq[1][1], analytical_freq[1], 0.1);
  BOOST_CHECK_CLOSE(t[1], 2*dt, 0.1);
  BOOST_REQUIRE_EQUAL(instant_freq[2].size(), 2);
  BOOST_CHECK_CLOSE(instant_freq[2][0], analytical_freq[0], 0.1);
  BOOST_CHECK_CLOSE(instant_freq[2][1], analytical_freq[1], 0.1);
  BOOST_CHECK_CLOSE(t[2], 3*dt, 0.1);
}


BOOST_FIXTURE_TEST_CASE(instantaneous_frequency_mean_field_and_phase_network,
    MKuramoto) {
  system->Integrate(dt, n, InstantaneousFrequencyMeanFieldAndPhaseObserver<
      MKuramotoSakaguchiSystem, state_type, std::vector<state_type>>((*system), 
        dt, 1, instant_freq, instant_freq_mean_field, t));
  BOOST_REQUIRE_EQUAL(instant_freq_mean_field.size(), n+1);
  BOOST_REQUIRE_EQUAL(instant_freq.size(), n+1);
  BOOST_REQUIRE_EQUAL(t.size(), n+1);
  BOOST_REQUIRE_EQUAL(instant_freq_mean_field[0].size(), 2);
  BOOST_CHECK_CLOSE(instant_freq_mean_field[0][0], analytical_freq[0], 0.1);
  BOOST_CHECK_CLOSE(instant_freq_mean_field[0][1], analytical_freq[1], 0.1);
  BOOST_REQUIRE_EQUAL(instant_freq[0].size(), 4);
  BOOST_CHECK_CLOSE(instant_freq[0][0], frequency[0], 0.01);
  BOOST_CHECK_CLOSE(instant_freq[0][1], frequency[1], 0.01);
  BOOST_CHECK_CLOSE(instant_freq[0][2], frequency[2], 0.01);
  BOOST_CHECK_CLOSE(instant_freq[0][3], frequency[3], 0.01);
  BOOST_CHECK_CLOSE(t[0], dt, 0.1);
  BOOST_REQUIRE_EQUAL(instant_freq_mean_field[1].size(), 2);
  BOOST_CHECK_CLOSE(instant_freq_mean_field[1][0], analytical_freq[0], 0.1);
  BOOST_CHECK_CLOSE(instant_freq_mean_field[1][1], analytical_freq[1], 0.1);
  BOOST_REQUIRE_EQUAL(instant_freq[1].size(), 4);
  BOOST_CHECK_CLOSE(instant_freq[1][0], frequency[0], 0.01);
  BOOST_CHECK_CLOSE(instant_freq[1][1], frequency[1], 0.01);
  BOOST_CHECK_CLOSE(instant_freq[1][2], frequency[2], 0.01);
  BOOST_CHECK_CLOSE(instant_freq[1][3], frequency[3], 0.01);
  BOOST_CHECK_CLOSE(t[1], 2*dt, 0.1);
  BOOST_REQUIRE_EQUAL(instant_freq_mean_field[2].size(), 2);
  BOOST_CHECK_CLOSE(instant_freq_mean_field[2][0], analytical_freq[0], 0.1);
  BOOST_CHECK_CLOSE(instant_freq_mean_field[2][1], analytical_freq[1], 0.1);
  BOOST_REQUIRE_EQUAL(instant_freq[2].size(), 4);
  BOOST_CHECK_CLOSE(instant_freq[2][0], frequency[0], 0.01);
  BOOST_CHECK_CLOSE(instant_freq[2][1], frequency[1], 0.01);
  BOOST_CHECK_CLOSE(instant_freq[2][2], frequency[2], 0.01);
  BOOST_CHECK_CLOSE(instant_freq[2][3], frequency[3], 0.01);
  BOOST_CHECK_CLOSE(t[2], 3*dt, 0.1);
}
