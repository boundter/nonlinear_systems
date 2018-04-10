#ifndef __AVERAGE_FREQUENCY_OBSERVER__
#define __AVERAGE_FREQUENCY_OBSERVER__

#include <vector>
#include <cmath>

namespace nonlinear_systems {
template <typename state_type = std::vector<double> >
class AverageFrequencyPhaseObserver {
  public:
    std::vector<double>& _average_frequency;

    AverageFrequencyPhaseObserver(const state_type& one_step_before,
        const state_type& two_steps_before, double dt,
        std::vector<double>& average_frequency)
      : _average_frequency(average_frequency) {
        _one_step_before = one_step_before;
        _two_steps_before = two_steps_before;
        _dt = dt;
        step_number = 0;
      }


    // TODO: Check length
    // TODO: Refactor moving average
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


template <typename state_type = std::vector<double> >
class AverageFrequencyObserver {
  public:
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
