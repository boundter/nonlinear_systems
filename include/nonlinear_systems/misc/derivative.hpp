#ifndef __DERIVATIVE__
#define __DERIVATIVE__

#include <cmath>
#include <vector>

namespace nonlinear_systems {
/*!
 *  \brief Calculates the derivative using the two point method.
 *
 *  The two point method calculates the derivative between points with
 *  \f[ f'(t) = \frac{f(t+h) - f(t-h)}{2h}, \f]
 *  where \f$ h \f$ is the timestep. This has an accuracy of \f$O(h^2)\f$.
 *  This function takes equidistant timesteps.
 *
 *  @param previous_state state before the current point in time.
 *  @param next_state state after the current point in time.
 *  @param timestep the timestep between states.
 */
template<typename state_type>
state_type TwoPointDerivative(const state_type& previous_state, 
    const state_type& next_state, double timestep) {
  // TODO: Check same dimension
  state_type derivative;
  double double_timestep = 2.*timestep;
  for (size_t i = 0; i < previous_state.size(); ++i) {
    double phase_difference = next_state[i] - previous_state[i];
    derivative.push_back(phase_difference/double_timestep);
  }
  return derivative;
}


/*!
 *  \brief Calculates the derivative using the two point method with modulo.
 *
 *  The two point method calculates the derivative between points with
 *  \f[ f'(t) = \frac{f(t+h) - f(t-h)}{2h}, \f]
 *  where \f$ h \f$ is the timestep. This has an accuracy of \f$O(h^2)\f$.
 *  This function takes equidistant timesteps and unwraps derivatives that take
 *  a bigger step than limit_step_size with the factor modulo, by adding it, if
 *  the previous state was positive and subtractinmg it, if it was negative. 
 *
 *  @param previous_state state before the current point in time.
 *  @param next_state state after the current point in time.
 *  @param timestep the timestep between states.
 *  @param modulo factor used to unwrap the state.
 *  @param limit_step_size upper limit for the magnitude of the steps, before
 *  they are unwrapped.
 */
template<typename state_type>
state_type TwoPointDerivative(const state_type& previous_state, 
    const state_type& next_state, double timestep, double modulo, 
    double limit_step_size = 1) {
  // TODO: Check same dimension
  state_type derivative;
  double double_timestep = 2.*timestep;
  for (size_t i = 0; i < previous_state.size(); ++i) {
    double phase_difference = next_state[i] - previous_state[i];
    if (fabs(phase_difference) > limit_step_size) {
      phase_difference += std::copysign(modulo, previous_state[i]);
    }
    derivative.push_back(phase_difference/double_timestep);
  }
  return derivative;
}
} // nonlinear_systems
#endif
