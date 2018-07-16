#define BOOST_TEST_MODULE ReducedMKuramotoSakaguchiWatanabeStrogatzSystem
#include <boost/test/included/unit_test.hpp>
#include <nonlinear_systems/systems/reduced_m_kuramoto_sakaguchi_watanabe_strogatz_system.hpp>

typedef std::vector<double> state_type;
typedef std::vector<state_type> network_type;

using namespace nonlinear_systems;

BOOST_AUTO_TEST_CASE(conversion_to_constants) {
  network_type phases = {{0.1, -0.2*M_PI, -0.2, 0.5},
                         {1.3, 1.5, M_PI, M_PI/2.}};
  ReducedMKuramotoSakaguchiWatanabeStrogatzSystem system(phases);
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
