#define BOOST_TEST_MODULE MKuramotoSakaguchiODE
#include <boost/test/included/unit_test.hpp>
#include <vector>
#include <nonlinear_systems/odes/m_kuramoto_sakaguchi_ode.hpp>

typedef std::vector<double> state_type;
typedef std::vector<state_type> network_type;
typedef std::vector<unsigned int> node_size_type;

using namespace nonlinear_systems;


struct F {
  F() {
    node_indices = {0, 2, 4};
    frequency = {0., 0.5, 1., 2.};
    state_type coupling_1 = {1., -1.};
    state_type coupling_2 = {0.5, 0.3};
    coupling = {coupling_1, coupling_2};
    state_type phase_shift_1 = {0.25, 0.15};
    state_type phase_shift_2 = {-0.15, 0.05};
    phase_shift = {phase_shift_1, phase_shift_2};
  }

  ~ F() {;}
  
  state_type frequency;
  network_type coupling;
  network_type phase_shift;
  node_size_type node_indices;
};


BOOST_FIXTURE_TEST_CASE(constructor, F) {
  // test for correct initalization
  MKuramotoSakaguchiODE ode(frequency, coupling, phase_shift, node_indices);
  state_type _frequency = ode._frequency;
  BOOST_REQUIRE_EQUAL(_frequency.size(), frequency.size());
  for (size_t i = 0; i < _frequency.size(); ++i) {
    BOOST_CHECK_CLOSE(_frequency[i], frequency[i], 0.01);
  }
  network_type _coupling = ode._coupling;
  BOOST_REQUIRE_EQUAL(_coupling.size(), coupling.size());
  for (size_t i = 0; i < _coupling.size(); ++i) {
    BOOST_REQUIRE_EQUAL(_coupling[i].size(), coupling[i].size());
    for (size_t j = 0; j < _coupling[i].size(); ++j) {
      BOOST_CHECK_CLOSE(_coupling[i][j], coupling[i][j], 0.01);
    }
  }
  network_type _phase_shift = ode._phase_shift;
  BOOST_REQUIRE_EQUAL(_phase_shift.size(), phase_shift.size());
  for (size_t i = 0; i < _phase_shift.size(); ++i) {
    BOOST_REQUIRE_EQUAL(_phase_shift[i].size(), phase_shift[i].size());
    for (size_t j = 0; j < _phase_shift[i].size(); ++j) {
      BOOST_CHECK_CLOSE(_phase_shift[i][j], phase_shift[i][j], 0.01);
    }
  }
  node_size_type _node_indices = ode._node_indices;
  BOOST_REQUIRE_EQUAL(_node_indices.size(), node_indices.size());
  for (size_t i = 0; i < _node_indices.size(); ++i) {
    BOOST_CHECK_EQUAL(_node_indices[i], node_indices[i]);
  }

  // test with different length
}

/*
BOOST_AUTO_TEST_CASE(mean_field) {

}


BOOST_AUTO_TEST_CASE(integration) {

}
*/
