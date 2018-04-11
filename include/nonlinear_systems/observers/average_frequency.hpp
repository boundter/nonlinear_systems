#ifndef __AVERAGE_FREQUENCY_OBSERVER__
#define __AVERAGE_FREQUENCY_OBSERVER__

#include <cmath> // copysign, atan2, fabs
#include <vector>

namespace nonlinear_systems {

/*!
 * This observer measures the average frequency of phase oscillators. It
 * calculates this for all oscillators using numerical differentiation.
 */
template <typename state_type = std::vector<double> >
class AverageFrequencyPhaseObserver {
  public:
    std::vector<double>& _average_frequency;

    
    /*!
     * @param one_step_before the phases one timestep before the integration
     * @param two_steps_before the phases two timesteps before the integration
     * @param dt the timestep of the integration
     * @param average_frequency a vector of the length of the phases in which
     * the frequencies will be saved
     */
    // TODO: Check length
    AverageFrequencyPhaseObserver(const state_type& one_step_before,
        const state_type& two_steps_before, double dt,
        std::vector<double>& average_frequency)
      : _average_frequency(average_frequency) {
        _one_step_before = one_step_before;
        _two_steps_before = two_steps_before;
        _dt = dt;
        step_number = 0;
      }


    // TODO: Refactor moving average
    // TODO: Refactor numerical differentiation
    // the differentiation uses the middlepoint method; the numerical derivative
    // x' at the timepoint n with timestep dt is:
    // x'[n] = (x[n+1] - x[n-1])/(2*dt)
    // the phases are temporarily unwrapped at the crossing from -pi to pi or pi
    // to -pi by subtracting/adding 2pi depending on the state before the
    // crossing
    void operator()(const state_type& x, double t) {
      for (size_t i = 0; i < _average_frequency.size(); ++i) {
        double phase_difference = x[i] - _two_steps_before[i];
        if (fabs(phase_difference) > 1.) {
          phase_difference += std::copysign(2*M_PI, _two_steps_before[i]);
        }
        double instant_frequency = phase_difference/(2.*_dt);
        _average_frequency[i] = 
          (instant_frequency + static_cast<double>(step_number)*_average_frequency[i])
          /(static_cast<double>(step_number) + 1);
      }
      _two_steps_before = _one_step_before;
      _one_step_before = x;
      step_number += 1;
    }
  
  
  protected:
    state_type _one_step_before;
    state_type _two_steps_before;
    double _dt;
    unsigned int step_number;
};


/*!
 *  This observer measures the average frequency of a limit cycle. For now it
 *  only works in two dimensions for states of the form [x1, y1, x2, y2, ...]
 *  and so on. The phase used is the normal polar coordinate atan2(y,x)
 */
template <typename state_type = std::vector<double> >
class AverageFrequencyObserver {
  public:

    /*!
     * @param one_step_before the phases one timestep before the integration
     * @param two_steps_before the phases two timesteps before the integration
     * @param dt the timestep of the integration
     * @param average_frequency a vector of the length of the phases in which
     * the frequencies will be saved
     */
    // TODO: Let user choose the coordinates to use
    AverageFrequencyObserver(const state_type& one_step_before,
        const state_type& two_steps_before, double dt,
        std::vector<double>& average_frequency) {
        std::vector<double> phases_one_before = CalculatePhases(one_step_before);
        std::vector<double> phases_two_before = CalculatePhases(two_steps_before);
        _phase_observer = 
          new AverageFrequencyPhaseObserver<std::vector<double> >(
              phases_one_before, phases_two_before, dt, average_frequency);
      }


    void operator()(const state_type& x, double t) {
      std::vector<double> phases = CalculatePhases(x);
      _phase_observer->operator()(phases, t);
    }


  protected:
    AverageFrequencyPhaseObserver<std::vector<double> >* _phase_observer;


    std::vector<double> CalculatePhases(const state_type& state) {
      std::vector<double> phases(state.size()/2);
      for (size_t i = 0; i < state.size()/2; ++i) {
        phases[i] = atan2(state[2*i+1], state[2*i]);
      }
      return phases;
    }
};
}

#endif
