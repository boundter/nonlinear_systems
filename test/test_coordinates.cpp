#define BOOST_TEST_MODULE Coordinates
#include <boost/test/included/unit_test.hpp>
#include <nonlinear_systems/misc/coordinates.hpp>
#include <vector>

typedef std::vector<double> state_type;

using namespace nonlinear_systems;

BOOST_AUTO_TEST_CASE(test_spherical) {
  // 2-dimensions
  state_type x_2 = {2., 3.};
  state_type spherical_2 = {3.606, 0.9828};
  state_type spherical = CartesianToSpherical(x_2);
  BOOST_TEST(spherical.size() == x_2.size());
  for (size_t i = 0; i < spherical.size(); ++i) {
    BOOST_CHECK_CLOSE_FRACTION(spherical[i], spherical_2[i], 0.01);
  }
  /*
  state_type spherical_iterator = CartesianToSpherical(x_2.begin(), x_2.end());
  BOOST_TEST(spherical_iterator.size() == x_2.size());
  for (size_t i = 0; i < spherical_iterator.size(); ++i) {
    BOOST_CHECK_CLOSE_FRACTION(spherical_iterator[i], spherical_2[i], 0.01);
  }
  */
  
  // 3-dimensions
  state_type x_3 = {6., 2., -1.};
  state_type spherical_3 = {6.403, 0.3567, -0.4636};
  spherical = CartesianToSpherical(x_3);
  BOOST_TEST(spherical.size() == x_3.size());
  for (size_t i = 0; i < spherical.size(); ++i) {
    BOOST_CHECK_CLOSE_FRACTION(spherical[i], spherical_3[i], 0.01);
  }
}
