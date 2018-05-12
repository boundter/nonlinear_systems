#ifndef __INSTANTANEOUS_FREQUENCY_OBSERVER__
#define __INSTANTANEOUS_FREQUENCY_OBSERVER__

#include <cmath>
#include <vector>
#include <nonlinear_systems/misc/derivative.hpp>

namespace nonlinear_systems {
template<typename system_type, typename state_type = std::vector<double>>
/*!
 *  \brief Observes the instantaneous frequency of all phases of the oscillators 
 *  of the system.
 */
class InstantaneousFrequencyObserver {
  public:
    std::vector<state_type>& _instantaneous_frequency;
    std::vector<double>& _t;
    system_type& _system;

    /*!
     *  During the construction the system will be integrated twice to generate
     *  the necessary history for the two-point derivative. This will give the
     *  first time of observation as t+dt.
     *
     *  @param system the observed system.
     *  @param dt the timestep.
     *  @param dimension the dimension of the oscillators.
     *  @param instantaneous_frequency the observed frequency.
     *  @param t the observed time.
     */
    InstantaneousFrequencyObserver(system_type& system, double dt, 
        unsigned int dimension, std::vector<state_type>& instantaneous_frequency,
        std::vector<double>& t) 
    : _instantaneous_frequency(instantaneous_frequency), _t(t), _system(system){
      _dt = dt;
      _d = dimension;
      _modulo = 2*M_PI;
      _limit_step_size = 1.;
      _two_steps_before = GetPhases(system.GetPositionSpherical());
      system.Integrate(dt, 1);
      _one_step_before = GetPhases(system.GetPositionSpherical());
      _t_before = system.GetTime();
      system.Integrate(dt, 1);
    }


    void operator()(const state_type& x, double t) {
      state_type current_state = _system.GetPositionSpherical();
      state_type phases = GetPhases(current_state);
      _instantaneous_frequency.push_back(TwoPointDerivative(_two_steps_before, 
            phases, _dt, _modulo, _limit_step_size));
      _two_steps_before = _one_step_before;
      _one_step_before = phases;
      _t.push_back(_t_before);
      _t_before = t;
    }

  protected:
    double _t_before;
    state_type _one_step_before;
    state_type _two_steps_before;
    double _dt;
    unsigned int _d;
    double _modulo;
    double _limit_step_size;


    // filters the order parameter from the state
    state_type GetPhases(const state_type& x) {
      state_type phases;
      if (_d == 1) {
        phases = x;
      }
      else {
        // TODO: check for correct size
        unsigned int N = x.size()/_d;
        for (unsigned int i = 0; i < N; ++i) {
          for (unsigned int j = i*_d+1; j < (i+1)*_d; ++j) {
            phases.push_back(x[j]);
          }
        }
      }
      return phases;
    }
};
} // nonlinear_systems
#endif
