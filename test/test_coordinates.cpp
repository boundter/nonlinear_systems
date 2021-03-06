#define BOOST_TEST_MODULE Coordinates
#include <boost/test/included/unit_test.hpp>
#include <cmath>
#include <vector>
#include <nonlinear_systems/misc/coordinates.hpp>

typedef std::vector<double> state_type;

using namespace nonlinear_systems;


BOOST_AUTO_TEST_CASE(CartesianToSpherical2D) {
  state_type x = {2., 3.};
  state_type analytical = {3.606, 0.9828};
  state_type spherical = CartesianToSpherical(x);
  BOOST_REQUIRE_EQUAL(spherical.size(), x.size());
  for (size_t i = 0; i < spherical.size(); ++i) {
    BOOST_CHECK_CLOSE(spherical[i], analytical[i], 0.1);
  }
}


BOOST_AUTO_TEST_CASE(CartesianToSpherical3D) {
  state_type x = {6., 2., -1.};
  state_type analytical = {6.403, 0.3567, -0.4636};
  state_type spherical = CartesianToSpherical(x);
  BOOST_REQUIRE_EQUAL(spherical.size(), x.size());
  for (size_t i = 0; i < spherical.size(); ++i) {
    BOOST_CHECK_CLOSE(spherical[i], analytical[i], 0.1);
  }
}


BOOST_AUTO_TEST_CASE(cartesian_to_spherical_interators) {
  // states of different dimension (2, 3); (6, 2, -1) 
  state_type x = {2, 3, 6, 2, -1};
  state_type analytical_2 = {3.606, 0.9828};
  state_type analytical_3 = {6.403, 0.3567, -0.4636};
  state_type spherical_2 = CartesianToSpherical<state_type>(x.begin(), x.begin()+2);
  BOOST_REQUIRE_EQUAL(spherical_2.size(), analytical_2.size());
  BOOST_CHECK_CLOSE(spherical_2[0], analytical_2[0], 0.1);
  BOOST_CHECK_CLOSE(spherical_2[1], analytical_2[1], 0.1);
  state_type spherical_3 = CartesianToSpherical<state_type>(x.begin()+2, x.end());
  BOOST_REQUIRE_EQUAL(spherical_3.size(), analytical_3.size());
  BOOST_CHECK_CLOSE(spherical_3[0], analytical_3[0], 0.1);
  BOOST_CHECK_CLOSE(spherical_3[1], analytical_3[1], 0.1);
  BOOST_CHECK_CLOSE(spherical_3[2], analytical_3[2], 0.1);
  
}
