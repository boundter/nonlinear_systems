#define BOOST_TEST_MODULE ReducedMKuramotoSakaguchiWatanabeStrogatzODE
#include <boost/test/included/unit_test.hpp>
#include <cmath>
#include <vector>
#include <nonlinear_systems/odes/reduced_m_kuramoto_sakaguchi_watanabe_strogatz_ode.hpp>

typedef std::vector<double> state_type;
typedef std::vector<state_type> network_type;
typedef std::vector<unsigned int> node_size_type;

using namespace nonlinear_systems;


BOOST_AUTO_TEST_CASE(mean_field) {
  node_size_type node_indices = {0, 3, 6};
  state_type phase_shift = {0.25, 0.3};
  state_type coupling = {2., 3.};
  network_type constants = {{-M_PI/2., M_PI/2.}, {-2.*M_PI/3., 0., 2.*M_PI/3.}};
  state_type state = {0.5, M_PI/2., 0., 0.3, M_PI/4., M_PI/2.};
  state_type omega = {0., 0.};
  unsigned int N = 5;

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
