#define BOOST_TEST_MODULE HarmonicOscillatorODE
#include <boost/test/included/unit_test.hpp>
#include <vector>
#include <nonlinear_systems/odes/harmonic_oscillator_ode.hpp>

typedef std::vector<double> state_type;

using namespace nonlinear_systems;


struct F {
  F() {
    params[0] = 2.;
    ode = std::unique_ptr<HarmonicOscillatorODE>(
        new HarmonicOscillatorODE(params));
  }
  
  ~F(){;}

  double* params;
  std::unique_ptr<HarmonicOscillatorODE> ode;
};


BOOST_FIXTURE_TEST_CASE(constructor, F) {
  // check if omega was correctly set
  BOOST_CHECK_CLOSE(ode->_omega, params[0], 0.01);
}


BOOST_FIXTURE_TEST_CASE(integration, F) {
  // check if the ode is correct
  state_type x= {2., 1.3};
  state_type dx(2);
  ode->operator()(x, dx, 0.);
  state_type analytic = {1.3, -8.};

  BOOST_REQUIRE_EQUAL(dx.size(), analytic.size());
  BOOST_CHECK_CLOSE(dx[0], analytic[0], 0.01);
  BOOST_CHECK_CLOSE(dx[1], analytic[1], 0.01);
}
