#define BOOST_TEST_MODULE ReducedMKuramotoSakaguchiWatanabeStrogatzSystem
#include <boost/test/included/unit_test.hpp>
#include <nonlinear_systems/systems/reduced_m_kuramoto_sakaguchi_watanabe_strogatz_system.hpp>

typedef std::vector<double> state_type;
typedef std::vector<state_type> network_type;

using namespace nonlinear_systems;

BOOST_AUTO_TEST_CASE(unknown_conversion) {
  network_type phases = {{0.1, -0.2*M_PI, -0.2, 0.5},
                         {1.3, 1.5, M_PI, M_PI/2.}};
  BOOST_REQUIRE_THROW(
  ReducedMKuramotoSakaguchiWatanabeStrogatzSystem(phases, "foo"), 
  std::invalid_argument);
  BOOST_REQUIRE_NO_THROW(
  ReducedMKuramotoSakaguchiWatanabeStrogatzSystem(phases, "identity"));
  BOOST_REQUIRE_NO_THROW(
  ReducedMKuramotoSakaguchiWatanabeStrogatzSystem(phases, "splay"));
}


BOOST_AUTO_TEST_CASE(identity_conversion) {
  network_type phases = {{0.1, -0.2*M_PI, -0.2, 0.5},
                         {1.3, 1.5, M_PI, M_PI/2.}};
  ReducedMKuramotoSakaguchiWatanabeStrogatzSystem system(phases, "identity");
  state_type state = system.GetPosition();
  BOOST_REQUIRE_EQUAL(state.size(), 6);
  BOOST_CHECK_SMALL(state[0], 0.01);
  BOOST_CHECK_SMALL(state[1], 0.01);
  BOOST_CHECK_SMALL(state[2], 0.01);
  BOOST_CHECK_SMALL(state[3], 0.01);
  BOOST_CHECK_SMALL(state[4], 0.01);
  BOOST_CHECK_SMALL(state[5], 0.01);
}


BOOST_AUTO_TEST_CASE(back_conversion_identity) {
  network_type phases = {{0.1, -0.2*M_PI, -0.2, 0.5},
                         {1.3, 1.5, M_PI, M_PI/2.}};
  ReducedMKuramotoSakaguchiWatanabeStrogatzSystem system(phases, "identity");
  network_type new_phases = system.CalculatePhases();
  BOOST_REQUIRE_EQUAL(new_phases.size(), phases.size());
  BOOST_REQUIRE_EQUAL(new_phases[0].size(), phases[0].size());
  BOOST_CHECK_CLOSE(new_phases[0][0], phases[0][0], 0.1);
  BOOST_CHECK_CLOSE(new_phases[0][1], phases[0][1], 0.1);
  BOOST_CHECK_CLOSE(new_phases[0][2], phases[0][2], 0.1);
  BOOST_CHECK_CLOSE(new_phases[0][3], phases[0][3], 0.1);
  BOOST_REQUIRE_EQUAL(new_phases[1].size(), phases[1].size());
  BOOST_CHECK_CLOSE(new_phases[1][0], phases[1][0], 0.1);
  BOOST_CHECK_CLOSE(new_phases[1][1], phases[1][1], 0.1);
  BOOST_CHECK_CLOSE(new_phases[1][2], phases[1][2], 0.1);
  BOOST_CHECK_CLOSE(new_phases[1][3], phases[1][3], 0.1);
}

BOOST_AUTO_TEST_CASE(splay_conversion) {
  network_type phases = {{0.1, -0.2*M_PI, -0.2, 0.5},
                         {1.3, 1.5, M_PI, M_PI/2.}};
  ReducedMKuramotoSakaguchiWatanabeStrogatzSystem system(phases, "splay");
  state_type state = system.GetPosition();
  BOOST_REQUIRE_EQUAL(state.size(), 6);
  BOOST_CHECK_CLOSE(state[0], 0.74295, 0.1);
  BOOST_CHECK_CLOSE(state[1], 0.02973, 0.1);
  BOOST_CHECK_CLOSE(state[2], -0.04451, 0.1);
  BOOST_CHECK_CLOSE(state[3], 0.89072, 0.1);
  BOOST_CHECK_CLOSE(state[4], -0.45517, 0.1);
  BOOST_CHECK_CLOSE(state[5], 1.50626, 0.1);
}


BOOST_AUTO_TEST_CASE(back_conversion_splay) {
  network_type phases = {{0.1, -0.2*M_PI, -0.2, 0.5},
                         {1.3, 1.5, M_PI, M_PI/2.}};
  ReducedMKuramotoSakaguchiWatanabeStrogatzSystem system(phases, "splay");
  network_type new_phases = system.CalculatePhases();
  BOOST_REQUIRE_EQUAL(new_phases.size(), phases.size());
  BOOST_REQUIRE_EQUAL(new_phases[0].size(), phases[0].size());
  BOOST_CHECK_CLOSE(new_phases[0][0], phases[0][0], 0.1);
  BOOST_CHECK_CLOSE(new_phases[0][1], phases[0][1], 0.1);
  BOOST_CHECK_CLOSE(new_phases[0][2], phases[0][2], 0.1);
  BOOST_CHECK_CLOSE(new_phases[0][3], phases[0][3], 0.1);
  BOOST_REQUIRE_EQUAL(new_phases[1].size(), phases[1].size());
  BOOST_CHECK_CLOSE(new_phases[1][0], phases[1][0], 0.1);
  BOOST_CHECK_CLOSE(new_phases[1][1], phases[1][1], 0.1);
  BOOST_CHECK_CLOSE(new_phases[1][2], phases[1][2], 0.1);
  BOOST_CHECK_CLOSE(new_phases[1][3], phases[1][3], 0.1);
}
