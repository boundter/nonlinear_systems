#define BOOST_TEST_MODULE Observers
#include <boost/test/included/unit_test.hpp>
#include <cmath>
#include <vector>
#include <boost/numeric/odeint/stepper/runge_kutta4.hpp>
#include <nonlinear_systems/systems/generic_system.hpp>
#include <nonlinear_systems/odes/kuramoto_sakaguchi_ode.hpp>
#include <nonlinear_systems/systems/m_kuramoto_sakaguchi_system.hpp>
#include <nonlinear_systems/observers/average_frequency_observer.hpp>

typedef std::vector<double> state_type;
typedef std::vector<state_type> network_type;
typedef boost::numeric::odeint::runge_kutta4<state_type> stepper_type;

using namespace nonlinear_systems;

// TODO: Refactor simple odes
class SimpleRotators {
  public:
    double omega;

    SimpleRotators(void* params) {
      omega = reinterpret_cast<double*>(params)[0];
    }

    void operator()(const state_type& x, state_type& dx, const double t) {
      dx[0] = omega;
      dx[1] = 2*omega;
    }
};


class TwoHarmonicOscillators {
  public:
    double omega;

    TwoHarmonicOscillators(void* params) {
      omega = reinterpret_cast<double*>(params)[0];
    }

    void operator()(const state_type& x, state_type& dx, const double t) {
      dx[0] = x[1];
      dx[1] = -omega*omega*x[0]; 
      dx[2] = x[3];
      dx[3] = -2*omega*2*omega*x[2];
    }
};


// TODO: Refactor in seperate file
class KuramotoSystem: public GenericSystem<KuramotoSakaguchiODE, state_type> {
    public:
      KuramotoSystem(unsigned int N_oscillator, double frequency, 
          double coupling, double phase_shift)
      : GenericSystem<KuramotoSakaguchiODE, state_type>(N_oscillator, 1) {
        this->_ode = std::unique_ptr<KuramotoSakaguchiODE>(new 
            KuramotoSakaguchiODE(N_oscillator, state_type(N_oscillator, frequency), 
              coupling, phase_shift));
      }
};


struct Kuramoto {
  Kuramoto() {
    dt = 0.01;
    N = 2;
    n_transient = 1e3;
    n_mean = 100;
    frequency = 1.;
    coupling = 1.;
    phase_shift = 0.;
    system = std::unique_ptr<KuramotoSystem>(new KuramotoSystem(N, frequency, 
          coupling, phase_shift));
    // measure system states for the calculation of the frequency
    state_type x = {0., 1.};
    system->SetPosition(x);
    system->Integrate(dt, n_transient);
    state_two_before = system->GetPosition();
    mean_field_two_before = system->CalculateMeanFieldSpherical();
    system->Integrate(dt, 1);
    state_one_before = system->GetPosition();
    mean_field_one_before = system->CalculateMeanFieldSpherical();
    system->Integrate(dt, 1);
  }

  ~Kuramoto() {;}

  double dt;
  unsigned int n_transient;
  unsigned int n_mean;
  unsigned int N;
  double frequency;
  double coupling;
  double phase_shift;
  std::unique_ptr<KuramotoSystem> system;
  state_type state_two_before, state_one_before, mean_field_two_before,
             mean_field_one_before;
};


struct MKuramoto {
  MKuramoto() {
  dt = 0.01;
  n_transient = 1e3;
  n_mean = 100;
  frequency = 1.; 
  eps = 0.;
  node_size = {4, 5};
  system = std::unique_ptr<MKuramotoSakaguchiSystem>(
      new MKuramotoSakaguchiSystem(frequency, eps, node_size));
  }

  double dt;
  unsigned int n_transient;
  unsigned int n_mean;
  double frequency;
  double eps;
  std::vector<unsigned int> node_size;
  std::unique_ptr<MKuramotoSakaguchiSystem> system;
};

// TODO: can be done with kuramoto
BOOST_AUTO_TEST_CASE(unwrapped_phase) {
  double params[] = {1.};
  GenericSystem<SimpleRotators> system(2, 1, params);

  double dt = 0.01;
  unsigned int n_steps = 1e4;
  state_type average_frequency(2);
  state_type two_steps_before = system.GetPosition();
  system.Integrate(dt, 1);
  state_type one_step_before = system.GetPosition();
  system.Integrate(dt, 1);
  system.Integrate(dt, n_steps, 
      AverageFrequencyPhaseObserver<state_type>(one_step_before, two_steps_before, 
        dt, average_frequency));

  BOOST_CHECK_CLOSE(params[0], average_frequency[0], 0.01);
  BOOST_CHECK_CLOSE(2*params[0], average_frequency[1], 0.01);
}


BOOST_AUTO_TEST_CASE(limit_cycle) {
  double params = 1.;
  GenericSystem<TwoHarmonicOscillators> system(2, 2, &params);

  double dt = 0.01;
  unsigned int n = 1e4;
  state_type average_frequency(2);
  state_type initial_condition {0., 1., 0., 1.};
  system.SetPosition(initial_condition);
  state_type two_steps_before = system.GetPosition();
  system.Integrate(dt, 1);
  state_type one_step_before = system.GetPosition();
  system.Integrate(dt, 1);
  system.Integrate(dt, n, 
      AverageFrequencyObserver<state_type>(one_step_before, two_steps_before, 
        dt, average_frequency));

  BOOST_CHECK_CLOSE_FRACTION(params, fabs(average_frequency[0]), 0.01);
  BOOST_CHECK_CLOSE_FRACTION(2*params, fabs(average_frequency[1]), 0.01);
}


BOOST_FIXTURE_TEST_CASE(mean_field_phase, Kuramoto) {
  state_type average_frequency(1);
  system->Integrate(dt, n_mean, 
      AverageFrequencyMeanFieldObserver<KuramotoSystem, state_type, state_type>(
        (*system), mean_field_one_before, mean_field_two_before, dt, average_frequency));
  BOOST_REQUIRE_EQUAL(average_frequency.size(), 1);
  BOOST_CHECK_CLOSE_FRACTION(average_frequency[0], 1., 0.01);
}


BOOST_FIXTURE_TEST_CASE(mean_field_phase_network, MKuramoto) {
  system->Integrate(dt, n_transient);
  network_type mean_field_two_before = system->CalculateMeanField();
  system->Integrate(dt, 1);
  network_type mean_field_one_before = system->CalculateMeanField();
  system->Integrate(dt, 1);
  state_type average_frequency_network;
  system->Integrate(dt, n_mean,
      AverageFrequencyMeanFieldObserver<MKuramotoSakaguchiSystem,
      state_type, network_type >((*system), mean_field_one_before,
        mean_field_two_before, dt, average_frequency_network));
  state_type result = {-0.1206, 0.5293}; // not quite sure if this is correct
  BOOST_REQUIRE_EQUAL(average_frequency_network.size(), 2);
  BOOST_CHECK_CLOSE_FRACTION(average_frequency_network[0], result[0], 0.01);
  BOOST_CHECK_CLOSE_FRACTION(average_frequency_network[1], result[1], 0.01);
}


BOOST_FIXTURE_TEST_CASE(mean_field_and_phase, Kuramoto) {
  state_type average_frequency;
  state_type average_frequency_mean_field;
  system->Integrate(dt, n_mean,
      AverageFrequencyMeanFieldAndPhaseObserver<KuramotoSystem,
      state_type, state_type>((*system), state_one_before, state_two_before,
        mean_field_one_before, mean_field_two_before, dt, average_frequency,
        average_frequency_mean_field));
  BOOST_REQUIRE_EQUAL(average_frequency.size(), 2);
  BOOST_REQUIRE_EQUAL(average_frequency_mean_field.size(), 1);
  BOOST_CHECK_CLOSE_FRACTION(average_frequency[0], 1., 0.1);
  BOOST_CHECK_CLOSE_FRACTION(average_frequency[1], 1., 0.1);
  BOOST_CHECK_CLOSE_FRACTION(average_frequency_mean_field[0], 1., 0.1);
}


BOOST_FIXTURE_TEST_CASE(mean_field_and_phase_integration_in_constructor, 
    Kuramoto) {
  state_type average_frequency;
  state_type average_frequency_mean_field;
  system->Integrate(dt, n_mean,
      AverageFrequencyMeanFieldAndPhaseObserver<KuramotoSystem,
      state_type, state_type>((*system), dt, average_frequency, 
        average_frequency_mean_field));
  BOOST_REQUIRE_EQUAL(average_frequency.size(), 2);
  BOOST_REQUIRE_EQUAL(average_frequency_mean_field.size(), 1);
  BOOST_CHECK_CLOSE_FRACTION(average_frequency[0], 1., 0.1);
  BOOST_CHECK_CLOSE_FRACTION(average_frequency[1], 1., 0.1);
  BOOST_CHECK_CLOSE_FRACTION(average_frequency_mean_field[0], 1., 0.1);
}


BOOST_FIXTURE_TEST_CASE(mean_field_and_phase_network_integration_in_constructor,
    MKuramoto) {
  state_type average_frequency_phase;
  state_type average_frequency_network;
  network_type mean_field_0 = system->CalculateMeanField();
  system->Integrate(dt, n_mean,
      AverageFrequencyMeanFieldAndPhaseObserver<MKuramotoSakaguchiSystem,
      state_type, network_type>((*system), dt, average_frequency_phase, 
        average_frequency_network));
  network_type mean_field_1 = system->CalculateMeanField();
  state_type result;
  result.push_back((mean_field_1[0][1] - mean_field_0[0][1])/(dt*n_mean));
  result.push_back((mean_field_1[1][1] - mean_field_0[1][1])/(dt*n_mean));
  BOOST_REQUIRE_EQUAL(average_frequency_network.size(), 2);
  BOOST_CHECK_CLOSE_FRACTION(average_frequency_network[0], result[0], 0.1);
  BOOST_CHECK_CLOSE_FRACTION(average_frequency_network[1], result[1], 0.1);
}
