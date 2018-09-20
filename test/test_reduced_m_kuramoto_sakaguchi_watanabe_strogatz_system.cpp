#define BOOST_TEST_MODULE ReducedMKuramotoSakaguchiWatanabeStrogatzSystem
#include <boost/test/included/unit_test.hpp>
#include <cmath>
#include <nonlinear_systems/systems/reduced_m_kuramoto_sakaguchi_watanabe_strogatz_system.hpp>
#include <nonlinear_systems/systems/m_kuramoto_sakaguchi_system.hpp>

typedef std::vector<double> state_type;
typedef std::vector<state_type> network_type;

using namespace nonlinear_systems;


BOOST_AUTO_TEST_CASE(conversion) {
  network_type phases = {{0.1, -0.2*M_PI, -0.2, 0.5},
                         {1.3, 1.5, M_PI, M_PI/2.}};
  state_type coupling = {1., 2.};
  state_type phase_shift = {1., 2.};
  state_type frequency = {1., 2.};
  ReducedMKuramotoSakaguchiWatanabeStrogatzSystem system(phases, frequency, 
      coupling, phase_shift);
  state_type state = system.GetPosition();
  BOOST_REQUIRE_EQUAL(state.size(), 6);
  BOOST_CHECK_CLOSE(state[0], 0.74295, 0.1);
  BOOST_CHECK_CLOSE(state[1], 0.02973, 0.1);
  BOOST_CHECK_CLOSE(state[2], -0.04451, 0.1);
  BOOST_CHECK_CLOSE(state[3], 0.89072, 0.1);
  BOOST_CHECK_CLOSE(state[4], -0.45517, 0.1);
  BOOST_CHECK_CLOSE(state[5], 1.50626, 0.1);
}


BOOST_AUTO_TEST_CASE(back_conversion) {
  network_type phases = {{0.1, -0.2*M_PI, -0.2, 0.5},
                         {1.3, 1.5, M_PI, M_PI/2.}};
  state_type coupling = {1., 2.};
  state_type phase_shift = {1., 2.};
  state_type frequency = {1., 2.};
  ReducedMKuramotoSakaguchiWatanabeStrogatzSystem system(phases, frequency, 
      coupling, phase_shift);
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


BOOST_AUTO_TEST_CASE(simple_constructor) {
  double frequency = 1.;
  double repulsive_excess = 2.;
  node_size_type node_size = {4, 6};
  unsigned int seed = 123456789;
  double cluster_width = 0.01;
  double cluster_distance = 0.25*M_PI;
  BOOST_REQUIRE_THROW(
      ReducedMKuramotoSakaguchiWatanabeStrogatzSystem(frequency, 
        repulsive_excess, node_size, "foo", seed, cluster_width,
        cluster_distance), std::invalid_argument);
  BOOST_REQUIRE_NO_THROW(
      ReducedMKuramotoSakaguchiWatanabeStrogatzSystem(frequency, 
        repulsive_excess, node_size, "clusters", seed, cluster_width,
        cluster_distance));
  BOOST_REQUIRE_NO_THROW(
      ReducedMKuramotoSakaguchiWatanabeStrogatzSystem(frequency, 
        repulsive_excess, node_size, "random", seed, cluster_width,
        cluster_distance));
  ReducedMKuramotoSakaguchiWatanabeStrogatzSystem random_system(frequency, 
    repulsive_excess, node_size, "random");
  network_type phases = random_system.CalculatePhases();
  BOOST_REQUIRE_EQUAL(phases.size(), node_size.size());
  BOOST_REQUIRE_EQUAL(phases[0].size(), node_size[0]);
  BOOST_REQUIRE_EQUAL(phases[1].size(), node_size[1]);
  ReducedMKuramotoSakaguchiWatanabeStrogatzSystem cluster_system(frequency, 
    repulsive_excess, node_size, "clusters");
  phases = cluster_system.CalculatePhases();
  BOOST_REQUIRE_EQUAL(phases.size(), node_size.size());
  BOOST_REQUIRE_EQUAL(phases[0].size(), node_size[0]);
  BOOST_CHECK_SMALL(phases[0][0] - phases[0][1], 2*cluster_width);
  BOOST_CHECK_SMALL(phases[0][0] - phases[0][2], 2*cluster_width);
  BOOST_CHECK_SMALL(phases[0][0] - phases[0][3], 2*cluster_width);
  BOOST_REQUIRE_EQUAL(phases[1].size(), node_size[1]);
  BOOST_CHECK_SMALL(phases[1][0] - phases[1][1], 2*cluster_width);
  BOOST_CHECK_SMALL(phases[1][0] - phases[1][2], 2*cluster_width);
  BOOST_CHECK_SMALL(phases[1][0] - phases[1][3], 2*cluster_width);
  BOOST_CHECK_SMALL(phases[1][0] - phases[1][4], 2*cluster_width);
  BOOST_CHECK_SMALL(phases[1][0] - phases[1][5], 2*cluster_width);
}


BOOST_AUTO_TEST_CASE(integration_splay_random) {
  double frequency = 0.5;
  double repulsive_excess = 0.1;
  node_size_type node_size = {4, 6};
  unsigned int seed = 123456789;
  double cluster_width = 0.1;
  double cluster_distance = 0.25*M_PI;
  
  ReducedMKuramotoSakaguchiWatanabeStrogatzSystem ws_system(frequency,
      repulsive_excess, node_size, "random", seed);
  MKuramotoSakaguchiSystem kuramoto_system(frequency, repulsive_excess, 
      node_size, seed);
  // Same initial conditions
  network_type ws_phases = ws_system.CalculatePhases();
  network_type kuramoto_phases = kuramoto_system.GetNodes();
  BOOST_REQUIRE_EQUAL(ws_phases.size(), kuramoto_phases.size());
  BOOST_REQUIRE_EQUAL(ws_phases[0].size(), kuramoto_phases[0].size());
  BOOST_CHECK_CLOSE(ws_phases[0][0], kuramoto_phases[0][0], 0.1);
  BOOST_CHECK_CLOSE(ws_phases[0][1], kuramoto_phases[0][1], 0.1);
  BOOST_CHECK_CLOSE(ws_phases[0][2], kuramoto_phases[0][2], 0.1);
  BOOST_CHECK_CLOSE(ws_phases[0][3], kuramoto_phases[0][3], 0.1);
  BOOST_REQUIRE_EQUAL(ws_phases[1].size(), kuramoto_phases[1].size());
  BOOST_CHECK_CLOSE(ws_phases[1][0], kuramoto_phases[1][0], 0.1);
  BOOST_CHECK_CLOSE(ws_phases[1][1], kuramoto_phases[1][1], 0.1);
  BOOST_CHECK_CLOSE(ws_phases[1][2], kuramoto_phases[1][2], 0.1);
  BOOST_CHECK_CLOSE(ws_phases[1][3], kuramoto_phases[1][3], 0.1);
  BOOST_CHECK_CLOSE(ws_phases[1][4], kuramoto_phases[1][4], 0.1);
  BOOST_CHECK_CLOSE(ws_phases[1][5], kuramoto_phases[1][5], 0.1);
  
  // Same final state
  double dt = 0.01;
  unsigned int number_steps = 1e3;
  ws_system.Integrate(dt, number_steps);
  kuramoto_system.Integrate(dt, number_steps);
  ws_phases = ws_system.CalculatePhases();
  kuramoto_phases = kuramoto_system.GetNodes();
  BOOST_REQUIRE_EQUAL(ws_phases.size(), kuramoto_phases.size());
  BOOST_REQUIRE_EQUAL(ws_phases[0].size(), kuramoto_phases[0].size());
  BOOST_CHECK_CLOSE(fmod(ws_phases[0][0], 2*M_PI), 
      fmod(kuramoto_phases[0][0] + 2*M_PI, 2*M_PI), 0.1);
  BOOST_CHECK_CLOSE(fmod(ws_phases[0][1], 2*M_PI), 
      fmod(kuramoto_phases[0][1] + 2*M_PI, 2*M_PI), 0.1);
  BOOST_CHECK_CLOSE(fmod(ws_phases[0][2], 2*M_PI), 
      fmod(kuramoto_phases[0][2] + 2*M_PI, 2*M_PI), 0.1);
  BOOST_CHECK_CLOSE(fmod(ws_phases[0][3], 2*M_PI), 
      fmod(kuramoto_phases[0][3] + 2*M_PI, 2*M_PI), 0.1);
  BOOST_REQUIRE_EQUAL(ws_phases[1].size(), kuramoto_phases[1].size());
  BOOST_CHECK_CLOSE(fmod(ws_phases[1][0] + 2*M_PI, 2*M_PI), 
      fmod(kuramoto_phases[1][0], 2*M_PI), 0.1);
  BOOST_CHECK_CLOSE(fmod(ws_phases[1][1] + 2*M_PI, 2*M_PI), 
      fmod(kuramoto_phases[1][1], 2*M_PI), 0.1);
  BOOST_CHECK_CLOSE(fmod(ws_phases[1][2] + 2*M_PI, 2*M_PI), 
      fmod(kuramoto_phases[1][2], 2*M_PI), 0.1);
  BOOST_CHECK_CLOSE(fmod(ws_phases[1][3] + 2*M_PI, 2*M_PI), 
      fmod(kuramoto_phases[1][3], 2*M_PI), 0.1);
  BOOST_CHECK_CLOSE(fmod(ws_phases[1][4] + 2*M_PI, 2*M_PI), 
      fmod(kuramoto_phases[1][4], 2*M_PI), 0.1);
  BOOST_CHECK_CLOSE(fmod(ws_phases[1][5] + 2*M_PI, 2*M_PI), 
      fmod(kuramoto_phases[1][5], 2*M_PI), 0.1);
}
