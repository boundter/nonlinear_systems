#define BOOST_TEST_MODULE Coordinates
#include <boost/test/included/unit_test.hpp>
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
