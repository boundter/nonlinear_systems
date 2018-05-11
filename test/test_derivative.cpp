#define BOOST_TEST_MODULE Derivative
#include <boost/test/included/unit_test.hpp>
#include <vector>
#include <nonlinear_systems/misc/derivative.hpp>

typedef std::vector<double> state_type;

using namespace nonlinear_systems;


struct F {
  F() {
    previous_state = {1., 2., 1.};
    next_state = {2., 3., 9.};
    dt = 0.25;
  }

  state_type previous_state;
  state_type next_state;
  double dt;
};

BOOST_FIXTURE_TEST_CASE(two_point_derivative, F) {
  state_type analytical = {2., 2., 16.};
  state_type derivative = TwoPointDerivative(previous_state, next_state, dt);
  BOOST_REQUIRE_EQUAL(derivative.size(), analytical.size());
  BOOST_CHECK_CLOSE(derivative[0], analytical[0], 0.01);
  BOOST_CHECK_CLOSE(derivative[1], analytical[1], 0.01);
  BOOST_CHECK_CLOSE(derivative[2], analytical[2], 0.01);
}


BOOST_FIXTURE_TEST_CASE(two_point_derivative_modulo, F) {
  double modulo = 6;
  next_state = {2., 3., 3.};
  double limit_step_size = 1.5;
  state_type analytical = {2., 2., 16.};
  state_type derivative = TwoPointDerivative(previous_state, next_state, dt, 
      modulo, limit_step_size);
  BOOST_REQUIRE_EQUAL(derivative.size(), analytical.size());
  BOOST_CHECK_CLOSE(derivative[0], analytical[0], 0.01);
  BOOST_CHECK_CLOSE(derivative[1], analytical[1], 0.01);
  BOOST_CHECK_CLOSE(derivative[2], analytical[2], 0.01);
}
