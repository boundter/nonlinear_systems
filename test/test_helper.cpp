#define BOOST_TEST_MODULE Helper
#include <boost/test/included/unit_test.hpp>
#include <vector>
#include <nonlinear_systems/misc/helper.hpp>

typedef std::vector<double> state_type;
typedef std::vector<state_type> network_type;

using namespace nonlinear_systems;

struct F {
  F() {
    state = {5.3, 2.1, 5.4, 2.1};
    network = {{1.2, 3.4, 4.5}, {6.3, 1.2, 7.8}};
  }

  ~F() {;}

  state_type state;
  network_type network;
  MeanFieldHelper<state_type> mean_field_helper_state;
  MeanFieldHelper<network_type> mean_field_helper_network;
};

BOOST_FIXTURE_TEST_CASE(phase, F) {
  // vector
  state_type phase_state = mean_field_helper_state.GetMeanFieldPhase(state);
  BOOST_REQUIRE_EQUAL(phase_state.size(), state.size() - 1);
  BOOST_CHECK_CLOSE(phase_state[0], state[1], 0.01);
  BOOST_CHECK_CLOSE(phase_state[1], state[2], 0.01);
  BOOST_CHECK_CLOSE(phase_state[2], state[3], 0.01);
  // matrix
  state_type phase_network = mean_field_helper_network.GetMeanFieldPhase(network);
  BOOST_REQUIRE_EQUAL(phase_network.size(), network[0].size()-1 + network[1].size()-1);
  BOOST_CHECK_CLOSE(phase_network[0], network[0][1], 0.01);
  BOOST_CHECK_CLOSE(phase_network[1], network[0][2], 0.01);
  BOOST_CHECK_CLOSE(phase_network[2], network[1][1], 0.01);
  BOOST_CHECK_CLOSE(phase_network[3], network[1][2], 0.01);
}


BOOST_FIXTURE_TEST_CASE(order, F) {
  // vector
  state_type order_state = mean_field_helper_state.GetOrderParameter(state);
  BOOST_REQUIRE_EQUAL(order_state.size(), 1);
  BOOST_CHECK_CLOSE(order_state[0], state[0], 0.01);
  // matrix
  state_type order_network = mean_field_helper_network.GetOrderParameter(network);
  BOOST_REQUIRE_EQUAL(order_network.size(), network.size());
  BOOST_CHECK_CLOSE(order_network[0], network[0][0], 0.01);
  BOOST_CHECK_CLOSE(order_network[1], network[1][0], 0.01);
}
