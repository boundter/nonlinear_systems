#define BOOST_TEST_MODULE KuramotoSakaguchiODE
#include <boost/test/included/unit_test.hpp>
#include <cmath>
#include <stdexcept>
#include <vector>
#include <nonlinear_systems/odes/kuramoto_sakaguchi_ode.hpp>

typedef std::vector<double> state_type;

using namespace nonlinear_systems;

struct F {
  F() {
    N = 3;
    frequency = {1., 1.2, 2.};
    coupling = 2.;
    phase_shift = 0.5;
    ode = std::unique_ptr<KuramotoSakaguchiODE>(new KuramotoSakaguchiODE(
          N, frequency, coupling, phase_shift));
  }

  ~F() {;}

  unsigned int N;
  state_type frequency;
  double coupling;
  double phase_shift;
  std::unique_ptr<KuramotoSakaguchiODE> ode;
};


BOOST_FIXTURE_TEST_CASE(constructor, F) {
  // check for correct length, too short and too long
  BOOST_REQUIRE_NO_THROW(KuramotoSakaguchiODE(N, frequency, coupling, phase_shift));
  BOOST_CHECK_THROW(KuramotoSakaguchiODE(N, {0., 0.}, coupling, phase_shift),
      std::length_error);
  BOOST_CHECK_THROW(KuramotoSakaguchiODE(N, {0., 0., 0., 0.}, coupling, 
      phase_shift), std::length_error);

  // check for coorect initalization
  BOOST_REQUIRE_EQUAL(ode->_N, N);
  BOOST_REQUIRE_EQUAL(ode->_frequency.size(), frequency.size());
  BOOST_CHECK_CLOSE(ode->_frequency[0], frequency[0], 0.01);
  BOOST_CHECK_CLOSE(ode->_frequency[1], frequency[1], 0.01);
  BOOST_CHECK_CLOSE(ode->_frequency[2], frequency[2], 0.01);
  BOOST_CHECK_CLOSE(ode->_coupling, coupling, 0.01);
  BOOST_CHECK_CLOSE(ode->_phase_shift, phase_shift, 0.01);
}


BOOST_FIXTURE_TEST_CASE(calculate_mean_field, F) {
  state_type x = {0., M_PI, M_PI/2.};
  state_type analytical = {0.3333, M_PI/2.};
  state_type mean_field = ode->CalculateMeanField(x);
  BOOST_REQUIRE_EQUAL(mean_field.size(), 2);
  BOOST_CHECK_CLOSE(mean_field[0], analytical[0], 0.1);
  BOOST_CHECK_CLOSE(mean_field[1], analytical[1], 0.1);
}


BOOST_FIXTURE_TEST_CASE(integration, F) {
  state_type x = {0., M_PI, M_PI/2.};
  state_type analytical = {1.585, 0.615, 2.319};
  state_type dx(N);
  ode->operator()(x, dx, 0.);
  BOOST_REQUIRE_EQUAL(dx.size(), analytical.size());
  BOOST_CHECK_CLOSE(dx[0], analytical[0], 0.1);
  BOOST_CHECK_CLOSE(dx[1], analytical[1], 0.1);
  BOOST_CHECK_CLOSE(dx[2], analytical[2], 0.1);
}
