#define BOOST_TEST_MODULE RayleighMeanFieldSystem
#include <boost/test/included/unit_test.hpp>
#include <nonlinear_systems/odes/rayleigh_mean_field_ode.hpp>
#include <nonlinear_systems/systems/rayleigh_mean_field_system.hpp>
#include <stdexcept>
#include <cmath>

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
  RayleighMeanFieldODEX ode = RayleighMeanFieldODEX(N, frequency, nonlinearity,
      coupling);

  // Test for wrong length of frequency
  BOOST_CHECK_THROW(ode = RayleighMeanFieldODEX(N+1, frequency, nonlinearity, 
        coupling), std::length_error);
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
  RayleighMeanFieldODEY ode = RayleighMeanFieldODEY(N, frequency, nonlinearity,
      coupling);

  // Test for wrong length of frequency
  BOOST_CHECK_THROW(ode = RayleighMeanFieldODEY(N+1, frequency, nonlinearity, 
        coupling), std::length_error);
}
