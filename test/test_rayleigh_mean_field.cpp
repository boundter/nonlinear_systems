#define BOOST_TEST_MODULE RayleighMeanFieldSystem
#include <boost/test/included/unit_test.hpp>
#include <stdexcept>
#include <cmath>
#include <random>
#include <nonlinear_systems/odes/rayleigh_mean_field_ode.hpp>
#include <nonlinear_systems/systems/rayleigh_mean_field_system.hpp>

typedef std::vector<double> state_type;

using namespace nonlinear_systems;

// TODO: Add Tests for the integration
BOOST_AUTO_TEST_CASE(test_ODEX_constructor) {
  unsigned int N = 10;
  double nonlinearity = 6.;
  double coupling = 0.134;
  
  state_type frequency;
  for (unsigned int i = 0; i < N; ++i) {
    frequency.push_back(1.);
  }
  
  // Test for correct initialization
  BOOST_CHECK_NO_THROW(RayleighMeanFieldODEX(N, frequency, nonlinearity, coupling));

  // Test for wrong length of frequency
  BOOST_CHECK_THROW(RayleighMeanFieldODEX(N+1, frequency, nonlinearity, coupling), 
      std::length_error);
}


BOOST_AUTO_TEST_CASE(test_ODEY_constructor) {
  unsigned int N = 10;
  double nonlinearity = 6.;
  double coupling = 0.134;
  
  state_type frequency;
  for (unsigned int i = 0; i < N; ++i) {
    frequency.push_back(1.);
  }
  
  // Test for correct initialization
  BOOST_CHECK_NO_THROW(RayleighMeanFieldODEY(N, frequency, nonlinearity, coupling));

  // Test for wrong length of frequency
  BOOST_CHECK_THROW(RayleighMeanFieldODEY(N+1, frequency, nonlinearity, coupling), 
      std::length_error);
}

BOOST_AUTO_TEST_CASE(test_system_constructor) {
  unsigned int N = 10;
  double nonlinearity = 6.;
  double coupling = 0.134;

  // Test default constructor
  BOOST_CHECK_NO_THROW(RayleighMeanFieldSystem(N, nonlinearity, coupling));
  
  // Test constructor for x-coupling
  BOOST_CHECK_NO_THROW(RayleighMeanFieldSystem(N, nonlinearity, coupling, "x"));
  
  // Test constructor for y-coupling
  BOOST_CHECK_NO_THROW(RayleighMeanFieldSystem(N, nonlinearity, coupling, "y"));

  // Test constructor for invalid coupling
  BOOST_CHECK_THROW(RayleighMeanFieldSystem(N, nonlinearity, coupling, "u"),
      std::invalid_argument);

  // Check initialization of frequency and position from mt19937_64
  unsigned long int seed = 123456789;
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
  for (unsigned int i = 0; i < N; ++i) {
    frequency.push_back(normal(rng));
  }
  RayleighMeanFieldSystem system = RayleighMeanFieldSystem(N, nonlinearity,
      coupling, "y", seed);
  BOOST_TEST(system.GetPosition() == position, boost::test_tools::per_element());
  BOOST_TEST(system.GetFrequency() == frequency, boost::test_tools::per_element());

}
