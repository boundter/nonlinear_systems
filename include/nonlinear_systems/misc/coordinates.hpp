#ifndef __COORDINATES__
#define __COORDINATES__

#include <cmath>
#include <vector>

namespace nonlinear_systems {
/*!
 * \brief Calculates the coordinates on a sphere of the same dimension as
 * the phase space. If the dimension is 1, the corrdinates will be wrapped
 * around the unit circle. The first coordinate is the radius and the later
 * ones are the phases. Careful: in 3-d this is not the same as spherical
 * coordinates with polar angle and azimuth!
 */
template<typename state_type>
state_type CartesianToSpherical(const state_type& cartesian) {
  unsigned int dimension = cartesian.size();
  state_type spherical;
  state_type sum_squared(dimension);
  for (size_t i = 0; i < dimension; ++i) {
    sum_squared[0] += cartesian[i]*cartesian[i];
  }
  for (size_t i = 1; i < dimension; ++i) {
    sum_squared[i] = sum_squared[i-1] - cartesian[i-1]*cartesian[i-1];
  }
  // radius or distance to the origin
  spherical.push_back(sqrt(sum_squared[0]));
  // phases
  for (size_t i = 0; i < dimension - 2; ++i) {
    spherical.push_back(acos(cartesian[i]/sqrt(sum_squared[i])));
  }
  // last phase is in the range -pi/2 to pi/2
  if (cartesian.back() < 0) {
    spherical.push_back(-acos(cartesian[dimension-2]
                              /sqrt(sum_squared[dimension-2])));
  }
  else {
    spherical.push_back(acos(cartesian[dimension-2]
                             /sqrt(sum_squared[dimension-2])));
  }
  return spherical; 
}
} // nonlinear_systems

#endif
