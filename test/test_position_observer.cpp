#define BOOST_TEST_MODULE PositionObserver
#include <boost/test/included/unit_test.hpp>
#include <vector>
#include <nonlinear_systems/observers/position_observer.hpp>

typedef std::vector<double> state_type;

using namespace nonlinear_systems;


BOOST_AUTO_TEST_CASE(observation) {
  state_type x_1 = {1., 2., 3.};
  state_type x_2 = {5., 6., 7.};
  std::vector<state_type> observed_position;
  std::vector<double> observed_time;
  PositionObserver<state_type> observer(observed_position, observed_time);
  observer(x_1, 1.);
  observer(x_2, 2.);
  BOOST_REQUIRE_EQUAL(observed_position.size(), 2);
  BOOST_REQUIRE_EQUAL(observed_position[0].size(), 3);
  BOOST_REQUIRE_EQUAL(observed_position[1].size(), 3);
  BOOST_REQUIRE_EQUAL(observed_time.size(), 2);
  BOOST_CHECK_CLOSE(observed_position[0][0], x_1[0], 0.01);
  BOOST_CHECK_CLOSE(observed_position[0][1], x_1[1], 0.01);
  BOOST_CHECK_CLOSE(observed_position[0][2], x_1[2], 0.01);
  BOOST_CHECK_CLOSE(observed_position[1][0], x_2[0], 0.01);
  BOOST_CHECK_CLOSE(observed_position[1][1], x_2[1], 0.01);
  BOOST_CHECK_CLOSE(observed_position[1][2], x_2[2], 0.01);
  BOOST_CHECK_CLOSE(observed_time[0], 1., 0.01);
  BOOST_CHECK_CLOSE(observed_time[1], 2., 0.01);
}
