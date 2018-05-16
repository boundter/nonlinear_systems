#ifndef __FREQUENCY_HELPER__
#define __FREQUENCY_HELPER__

#include <cmath>
#include <vector>
#include <nonlinear_systems/misc/derivative.hpp>

namespace nonlinear_systems {
template<typename state_type = std::vector<double>>
class FrequencyHelper {
  public:
    // TODO: Check size
    // TODO: Test
    FrequencyHelper(const state_type& one_step_before, 
        const state_type& two_steps_before, double dt, double modulo, 
        double limit_step_size) {
      _one_step_before = one_step_before;
      _two_steps_before = two_steps_before;
      _dt = dt;
      _modulo = modulo;
      _limit_step_size = limit_step_size;
    }


    state_type operator()(const state_type& phases) {
      state_type frequency = TwoPointDerivative(_two_steps_before, phases, _dt,
          _modulo, _limit_step_size);
      _two_steps_before = _one_step_before;
      _one_step_before = phases;
      return frequency;
    }


  protected:
    state_type _one_step_before;    
    state_type _two_steps_before;
    double _dt;
    double _modulo;
    double _limit_step_size;
};
} // nonlinear_systems
#endif
