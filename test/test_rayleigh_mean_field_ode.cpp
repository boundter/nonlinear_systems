#define BOOST_TEST_MODULE RayleighMeanFieldODE
#include <boost/test/included/unit_test.hpp>
#include <stdexcept>
#include <vector>
#include <nonlinear_systems/odes/rayleigh_mean_field_ode.hpp>

typedef std::vector<double> state_type;

using namespace nonlinear_systems;


struct F {
  F() {
    N = 3;
    nonlinearity = 6.;
    coupling = 0.134;
    frequency = state_type(N, 2.);
  }

  ~F() {;}

  unsigned int N;
  double nonlinearity;
  double coupling;
  state_type frequency;
};


BOOST_FIXTURE_TEST_CASE(constructor_ODEX, F) {
  // correct length of frequency
  BOOST_CHECK_NO_THROW(RayleighMeanFieldODEX(N, frequency, nonlinearity, 
        coupling));
  // false length of frequency
  BOOST_CHECK_THROW(RayleighMeanFieldODEX(N+1, frequency, nonlinearity, 
        coupling), std::length_error);
  BOOST_CHECK_THROW(RayleighMeanFieldODEX(N-1, frequency, nonlinearity, 
        coupling), std::length_error);

  // correct initialization
  RayleighMeanFieldODEX ode(N, frequency, nonlinearity, coupling);
  BOOST_CHECK_EQUAL(N, ode._N);
  BOOST_CHECK_CLOSE(ode._nonlinearity, nonlinearity, 0.01);
  BOOST_CHECK_CLOSE(ode._coupling, coupling, 0.01);
  state_type omega = ode._frequency;
  BOOST_REQUIRE_EQUAL(omega.size(), frequency.size());
  BOOST_CHECK_CLOSE(omega[0], frequency[0], 0.01);
  BOOST_CHECK_CLOSE(omega[1], frequency[1], 0.01);
  BOOST_CHECK_CLOSE(omega[2], frequency[2], 0.01);

}


BOOST_FIXTURE_TEST_CASE(integration_ODEX, F) {
  state_type x = {0.5, 0.8, 1.2, -1.5, 0.2, 0.};
  state_type analytical = {0.8, -0.2541, -1.5, 6.374, 0., -0.741};
  state_type derivative(2*N);
  RayleighMeanFieldODEX ode(N, frequency, nonlinearity, coupling);
  ode(x, derivative, 0.);
  BOOST_REQUIRE_EQUAL(derivative.size(), analytical.size());
  BOOST_CHECK_CLOSE(derivative[0], analytical[0], 1);
  BOOST_CHECK_CLOSE(derivative[1], analytical[1], 1);
  BOOST_CHECK_CLOSE(derivative[2], analytical[2], 1);
  BOOST_CHECK_CLOSE(derivative[3], analytical[3], 1);
  BOOST_CHECK_SMALL(derivative[4], 0.1);
  BOOST_CHECK_CLOSE(derivative[5], analytical[5], 1);
}


BOOST_FIXTURE_TEST_CASE(constructor_ODEY, F) {
  // correct length of frequency
  BOOST_CHECK_NO_THROW(RayleighMeanFieldODEY(N, frequency, nonlinearity, 
        coupling));
  // false length of frequency
  BOOST_CHECK_THROW(RayleighMeanFieldODEY(N+1, frequency, nonlinearity, 
        coupling), std::length_error);
  BOOST_CHECK_THROW(RayleighMeanFieldODEY(N-1, frequency, nonlinearity, 
        coupling), std::length_error);
  
  // correct initialization
  RayleighMeanFieldODEY ode(N, frequency, nonlinearity, coupling);
  BOOST_CHECK_EQUAL(N, ode._N);
  BOOST_CHECK_CLOSE(ode._nonlinearity, nonlinearity, 0.01);
  BOOST_CHECK_CLOSE(ode._coupling, coupling, 0.01);
  state_type omega = ode._frequency;
  BOOST_REQUIRE_EQUAL(omega.size(), frequency.size());
  BOOST_CHECK_CLOSE(omega[0], frequency[0], 0.01);
  BOOST_CHECK_CLOSE(omega[1], frequency[1], 0.01);
  BOOST_CHECK_CLOSE(omega[2], frequency[2], 0.01);
}


BOOST_FIXTURE_TEST_CASE(integration_ODEY, F) {
  state_type x = {0.5, 0.8, 1.2, -1.5, 0.2, 0.};
  state_type analytical = {0.8, -0.4105, -1.5, 6.6197, 0., -0.831};
  state_type derivative(2*N);
  RayleighMeanFieldODEY ode(N, frequency, nonlinearity, coupling);
  ode(x, derivative, 0.);
  BOOST_REQUIRE_EQUAL(derivative.size(), analytical.size());
  BOOST_CHECK_CLOSE(derivative[0], analytical[0], 1);
  BOOST_CHECK_CLOSE(derivative[1], analytical[1], 1);
  BOOST_CHECK_CLOSE(derivative[2], analytical[2], 1);
  BOOST_CHECK_CLOSE(derivative[3], analytical[3], 1);
  BOOST_CHECK_SMALL(derivative[4], 0.1);
  BOOST_CHECK_CLOSE(derivative[5], analytical[5], 1);
}
