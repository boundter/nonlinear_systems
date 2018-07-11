#define BOOST_TEST_MODULE ReducedMKuramotoSakaguchiWatanabeStrogatzODE
#include <boost/test/included/unit_test.hpp>
#include <cmath>
#include <vector>
#include <nonlinear_systems/odes/reduced_m_kuramoto_sakaguchi_watanabe_strogatz_ode.hpp>

typedef std::vector<double> state_type;
typedef std::vector<state_type> network_type;
typedef std::vector<unsigned int> node_size_type;

using namespace nonlinear_systems;

struct F {
  F() {
    node_indices = {0, 3, 6};
    phase_shift = {0.25, 0.3};
    coupling = {2., 3.};
    // actually we need more than 2 oscillators in a node, but as a dummy this
    // works here
    constants = {{-M_PI/2., M_PI/2.}, {-2.*M_PI/3., 0., 2.*M_PI/3.}};
    state = {0.5, M_PI/2., 0., 0.3, M_PI/4., M_PI/2.};
    omega = {0., 1.};
    N = 5;
  }

  node_size_type node_indices;
  state_type phase_shift, coupling, state, omega;
  network_type constants;
  unsigned int N;
};


BOOST_FIXTURE_TEST_CASE(mean_field, F) {
  ReducedMKuramotoSakaguchiWatanabeStrogatzODE ode(constants, coupling, 
      phase_shift, omega, N, node_indices);
  network_type mean_field = ode.CalculateMeanField(state);
  network_type analytical_mean_field = {{-0.0253, 0.4484}, {0.4484, 0.0253}};
  BOOST_REQUIRE_EQUAL(mean_field.size(), analytical_mean_field.size());
  BOOST_REQUIRE_EQUAL(mean_field[0].size(), analytical_mean_field[0].size());
  BOOST_CHECK_CLOSE(mean_field[0][0], analytical_mean_field[0][0], 1);
  BOOST_CHECK_CLOSE(mean_field[0][1], analytical_mean_field[0][1], 1);
  BOOST_REQUIRE_EQUAL(mean_field[1].size(), analytical_mean_field[1].size());
  BOOST_CHECK_CLOSE(mean_field[1][0], analytical_mean_field[1][0], 1);
  BOOST_CHECK_CLOSE(mean_field[1][1], analytical_mean_field[1][1], 1);
}


BOOST_FIXTURE_TEST_CASE(integration, F) {
  ReducedMKuramotoSakaguchiWatanabeStrogatzODE ode(constants, coupling,
      phase_shift, omega, N, node_indices);
  state_type dx(6);
  ode(state, dx, 0.);
  state_type analytical = {-0.0094875, 0.3363, 0.5605, 0.204, 0.0384, 1.04596};
  BOOST_REQUIRE_EQUAL(dx.size(), analytical.size());
  BOOST_CHECK_CLOSE(dx[0], analytical[0], 1);
  BOOST_CHECK_CLOSE(dx[1], analytical[1], 0.1);
  BOOST_CHECK_CLOSE(dx[2], analytical[2], 0.1);
  BOOST_CHECK_CLOSE(dx[3], analytical[3], 0.1);
  BOOST_CHECK_CLOSE(dx[4], analytical[4], 0.1);
  BOOST_CHECK_CLOSE(dx[5], analytical[5], 0.1);
}
