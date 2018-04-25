#define BOOST_TEST_MODULE Observers
#include <boost/test/included/unit_test.hpp>
#include <nonlinear_systems/systems/generic_system.hpp>
#include <nonlinear_systems/observers/average_frequency.hpp>
#include <nonlinear_systems/systems/m_kuramoto_sakaguchi_system.hpp>
#include <cmath>
#include <iostream>

typedef std::vector<double> state_type;
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


class HarmonicOscillator {
  public:
    double omega;

    HarmonicOscillator(void* params) {
      omega = reinterpret_cast<double*>(params)[0];
    }

    void operator()(const state_type& x, state_type& dx, const double t) {
      dx[0] = x[1];
      dx[1] = -omega*omega*x[0]; 
      dx[2] = x[3];
      dx[3] = -2*omega*2*omega*x[2];
    }
};


class Kuramoto {
  public:
    double omega;
    double K;

    Kuramoto(void* params) {
      omega = reinterpret_cast<double*>(params)[0];
      K = reinterpret_cast<double*>(params)[1];
    }


    void operator()(const state_type& x, state_type& dx, double t) {
      state_type mean_field = CalculateMeanField(x);
      for (unsigned int i = 0; i < x.size(); ++i) {
        dx[i] = omega + mean_field[0]*K*sin(mean_field[1] - x[i]);
      }
    }

    state_type CalculateMeanField(const state_type& x) {
      state_type mean_field(2);
      double re = 0., im = 0.;
      for (unsigned int i = 0; i < x.size(); ++i) {
        re += cos(x[i]);
        im += sin(x[i]);
      }
      re /= static_cast<double>(x.size());
      im /= static_cast<double>(x.size());
      mean_field[0] = sqrt(re*re + im*im);
      mean_field[1] = atan2(im, re);
      return mean_field;
    }
};

class KuramotoSystem: public GenericSystem<Kuramoto, state_type> {
    public:
      KuramotoSystem(void* params, unsigned int N_oscillators)
      : GenericSystem<Kuramoto, state_type>(N_oscillators, 1, params) {};

      state_type CalculateMeanField() {
        return this->_ode->CalculateMeanField(this->_x);
      }
};


BOOST_AUTO_TEST_CASE(test_average_frequency_phase_unwrapped) {
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

// TODO: Why the wrong direction?
BOOST_AUTO_TEST_CASE(test_average_frequency) {
  double params[] = {1.};  
  GenericSystem<HarmonicOscillator> system(2, 2, params);

  double dt = 0.01;
  unsigned int n_steps = 1e4;
  state_type average_frequency(2);
  state_type initial_condition {0., 1., 0., 1.};
  system.SetPosition(initial_condition);
  state_type two_steps_before = system.GetPosition();
  system.Integrate(dt, 1);
  state_type one_step_before = system.GetPosition();
  system.Integrate(dt, 1);
  system.Integrate(dt, n_steps, 
      AverageFrequencyObserver<state_type>(one_step_before, two_steps_before, 
        dt, average_frequency));


  BOOST_CHECK_CLOSE_FRACTION(params[0], fabs(average_frequency[0]), 0.01);
  BOOST_CHECK_CLOSE_FRACTION(2*params[0], fabs(average_frequency[1]), 0.01);
}


BOOST_AUTO_TEST_CASE(test_mean_field_average_frequency) {
  double params[] = {1., 1.};
  KuramotoSystem system(params, 2);
  
  double dt = 0.01;
  unsigned int N_transient = 1e3;
  unsigned int N_mean = 100;
  state_type x = {0., 1.};
  system.SetPosition(x);
  system.Integrate(dt, N_transient);
  state_type state_two_before = system.CalculateMeanField();
  system.Integrate(dt, 1);
  state_type state_one_before = system.CalculateMeanField();
  system.Integrate(dt, 1);
  state_type average_frequency(1);
  system.Integrate(dt, N_mean, 
      AverageFrequencyMeanFieldPhaseObserver<KuramotoSystem, state_type, state_type>(
        system, state_one_before, state_two_before, dt, average_frequency));
  BOOST_CHECK_CLOSE_FRACTION(average_frequency[0], 1., 0.01);

  // Network
  double frequency = 1., eps = 0.;
  std::vector<unsigned int> node_size = {4, 5};
  MKuramotoSakaguchiSystem mkur(frequency, eps, node_size);
  dt = 0.01;
  N_transient = 1e3;
  N_mean = 100;
  mkur.Integrate(dt, N_transient);
  std::vector<std::vector<double> > mean_field_two_before = mkur.CalculateMeanField();
  mkur.Integrate(dt, 1);
  std::vector<std::vector<double> > mean_field_one_before = mkur.CalculateMeanField();
  mkur.Integrate(dt, 1);
  std::vector<double> average_frequency_network;
  mkur.Integrate(dt, N_mean,
      AverageFrequencyMeanFieldPhaseObserver<MKuramotoSakaguchiSystem,
      std::vector<double>, std::vector<std::vector<double> >
      >(mkur, mean_field_one_before,
        mean_field_two_before, dt, average_frequency_network));
  std::vector<double> result = {-0.1206, 0.5293}; // not quite sure if this is correct
  BOOST_TEST(average_frequency_network.size() == 2);
  BOOST_CHECK_CLOSE_FRACTION(average_frequency_network[0], result[0], 0.01);
  BOOST_CHECK_CLOSE_FRACTION(average_frequency_network[1], result[1], 0.01);
}


BOOST_AUTO_TEST_CASE(test_mean_field_and_phase_average_frequency) {
  double params[] = {1., 1.};
  KuramotoSystem system(params, 2);

  double dt = 0.01;
  unsigned int N_transient = 1e3;
  unsigned int N_mean = 100;
  state_type x = {0., 1.};
  system.SetPosition(x);
  system.Integrate(dt, N_transient);
  state_type state_two_before = system.GetPosition();
  state_type mean_field_two_before = system.CalculateMeanField();
  system.Integrate(dt, 1);
  state_type state_one_before = system.GetPosition();
  state_type mean_field_one_before = system.CalculateMeanField();
  system.Integrate(dt, 1);
  state_type average_frequency;
  state_type average_frequency_mean_field;
  system.Integrate(dt, N_mean,
      AverageFrequencyMeanFieldPhaseAndPhaseObserver<KuramotoSystem,
      state_type, state_type>(system, state_one_before, state_two_before,
        mean_field_one_before, mean_field_two_before, dt, average_frequency,
        average_frequency_mean_field));
  BOOST_TEST(average_frequency.size() == 2);
  BOOST_TEST(average_frequency_mean_field.size() == 1);
  BOOST_CHECK_CLOSE_FRACTION(average_frequency[0], 1., 0.1);
  BOOST_CHECK_CLOSE_FRACTION(average_frequency[1], 1., 0.1);
  BOOST_CHECK_CLOSE_FRACTION(average_frequency_mean_field[0], 1., 0.1);

  // Integration in observer
  state_type frequency;
  state_type frequency_mean_field;
  system.Integrate(dt, N_mean,
      AverageFrequencyMeanFieldPhaseAndPhaseObserver<KuramotoSystem,
      state_type, state_type>(system, dt, frequency, frequency_mean_field));
  BOOST_TEST(frequency.size() == 2);
  BOOST_TEST(frequency_mean_field.size() == 1);
  BOOST_CHECK_CLOSE_FRACTION(frequency[0], 1., 0.1);
  BOOST_CHECK_CLOSE_FRACTION(frequency[1], 1., 0.1);
  BOOST_CHECK_CLOSE_FRACTION(frequency_mean_field[0], 1., 0.1);

}
