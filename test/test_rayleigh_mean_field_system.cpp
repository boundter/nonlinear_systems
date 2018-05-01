#define BOOST_TEST_MODULE RayleighMeanFieldSystem
#include <boost/test/included/unit_test.hpp>
#include <stdexcept>
#include <cmath>
#include <random>
#include <nonlinear_systems/systems/rayleigh_mean_field_system.hpp>

typedef std::vector<double> state_type;

using namespace nonlinear_systems;

struct F {
  F() {
    N = 3;
    nonlinearity = 6.;
    coupling = 0.134;
  }

  ~F() {;}

  unsigned int N;
  double nonlinearity;
  double coupling;
};

BOOST_FIXTURE_TEST_CASE(constructor, F) {
  // Test constructor for y-coupling
  BOOST_CHECK_NO_THROW(RayleighMeanFieldSystem<RayleighMeanFieldODEY>(N, 
        nonlinearity, coupling));
  
  // Test constructor for x-coupling
  BOOST_CHECK_NO_THROW(RayleighMeanFieldSystem<RayleighMeanFieldODEX>(N, 
        nonlinearity, coupling));
  
  // Check initialization of frequency and position from mt19937_64
  unsigned long int seed = 123456789;
  RayleighMeanFieldSystem<RayleighMeanFieldODEY> system(N, nonlinearity,
      coupling, seed);
  double min_value = -3.;
  double max_value = 3.;
  double mean = 1.;
  double stdev = 0.01;
  std::mt19937_64 rng(seed);
  std::uniform_real_distribution<double> uniform(min_value, max_value);
  std::normal_distribution<double> normal(mean, stdev);
  state_type position, frequency;
  for (unsigned int i = 0; i < 2*N; ++i) {
    position.push_back(uniform(rng));
  }
  state_type x = system.GetPosition();
  BOOST_REQUIRE_EQUAL(x.size(), position.size());
  BOOST_CHECK_CLOSE(x[0], position[0], 0.01);
  BOOST_CHECK_CLOSE(x[1], position[1], 0.01);
  BOOST_CHECK_CLOSE(x[2], position[2], 0.01);
  for (unsigned int i = 0; i < N; ++i) {
    frequency.push_back(normal(rng));
  }
  state_type omega = system.GetFrequency();
  BOOST_REQUIRE_EQUAL(omega.size(), frequency.size());
  BOOST_CHECK_CLOSE(omega[0], frequency[0], 0.01);
  BOOST_CHECK_CLOSE(omega[1], frequency[1], 0.01);
  BOOST_CHECK_CLOSE(omega[2], frequency[2], 0.01);
}


BOOST_FIXTURE_TEST_CASE(period, F) {
  RayleighMeanFieldSystem<> system(N, nonlinearity,
      coupling);
  double dt = 0.01;
  double N_trans = 1e4;
  double N_int = 1e4;
  std::vector<state_type> states;
  std::vector<double> times;
  system.Integrate(dt, N_trans);
  double T = system.CalculateMeanPeriod(7, dt);
  double T_exp = 12.856;
  BOOST_CHECK_CLOSE(T, T_exp, 1);
}


BOOST_FIXTURE_TEST_CASE(setting_frequency, F) {
  RayleighMeanFieldSystem<> system(N, nonlinearity, coupling);

  state_type frequency_too_long(N+1);
  state_type frequency_too_short(N-1);
  state_type frequency(N);

  BOOST_CHECK_THROW(system.SetFrequency(frequency_too_long), std::length_error);
  BOOST_CHECK_THROW(system.SetFrequency(frequency_too_short), std::length_error);
  BOOST_CHECK_NO_THROW(system.SetFrequency(frequency));
}


BOOST_AUTO_TEST_CASE(setting_splay_state) {
  unsigned int N = 1000;
  double nonlinearity = 6;
  double coupling = 0.134;
  RayleighMeanFieldSystem<> system(N, nonlinearity, coupling);

  system.SetSplayState();
  
  // The limit cycle of the Rayleigh oscillator is point symetric to the origin,
  // so the mean field will vanish for this state
  state_type mean_field = system.CalculateMeanField();
  double order_parameter = sqrt(mean_field[0]*mean_field[0] 
                                + mean_field[1]*mean_field[1]);
  BOOST_CHECK_SMALL(order_parameter, 0.01);
}
