#ifndef __AVERAGE_FREQUENCY_OBSERVER__
#define __AVERAGE_FREQUENCY_OBSERVER__

#include <cmath> // copysign, atan2, fabs
#include <vector>
#include <nonlinear_systems/observers/statistics_observer.hpp>
#include <nonlinear_systems/misc/helper.hpp>
#include <nonlinear_systems/misc/frequency_helper.hpp>

typedef std::vector<double> frequency_type;


namespace nonlinear_systems {
/*!
 * /brief This observer measures the average frequency of phase oscillators. It
 * calculates this for all oscillators using numerical differentiation.
 */
template <typename system_type, typename state_type = std::vector<double> >
class AverageFrequencyObserver {
  public:
    system_type& _system;
    
    /*!
     * @system the observed system.
     * @param one_step_before the phases from GetPositionSpherical one 
     * timestep before the integration
     * @param two_steps_before the phases from GetPositionSpherical two 
     * timesteps before the integration
     * @param dt the timestep of the integration
     * @param dimension the dimension of the ODE.
     * @param average_frequency a vector of the length of the phases in which
     * the frequencies will be saved
     */
    // TODO: Check length
    AverageFrequencyObserver(system_type& system, const frequency_type& one_step_before,
        const frequency_type& two_steps_before, double dt, unsigned int dimension,
        frequency_type& average_frequency)
   : _system(system) {
        _d = dimension;
        average_frequency.resize(one_step_before.size());
        _average_observer = std::shared_ptr<AverageObserver<frequency_type> >(
            new AverageObserver<frequency_type>(average_frequency));
        double modulo = 2*M_PI;
        double limit_step_size = 1.;
        _frequency_helper = std::shared_ptr<FrequencyHelper<frequency_type>>(
          new FrequencyHelper<frequency_type>(one_step_before, two_steps_before, 
            dt, modulo, limit_step_size));
      }


    void operator()(const state_type& x, double t) {
      frequency_type phases = this->GetPhases();
      frequency_type frequency = _frequency_helper->operator()(phases);
      _average_observer->operator()(frequency, t);
    }
  
  
  protected:
    std::shared_ptr<AverageObserver<frequency_type>> _average_observer;
    std::shared_ptr<FrequencyHelper<frequency_type>> _frequency_helper;
    unsigned int _d;
    
    // filters the order parameter from the state
    virtual frequency_type GetPhases() {
      state_type x = _system.GetPositionSpherical();
      frequency_type phases;
      if (_d == 1) {
        for (size_t i = 0; i < x.size(); ++i) {
          phases.push_back(x[i]);
        }
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


// TODO: Add tests for non-phase oscillators
/*!
 * \brief This observer measures the average frequency of the mean field of a system of
 * phase oscillators.  It can be used for general systems as well as networks,
 * but the timestep has to be constant.
 */
template<typename system_type, typename state_type = std::vector<double>,
  typename mean_field_type = state_type>
class AverageFrequencyMeanFieldObserver:
  public AverageFrequencyObserver<system_type, frequency_type> {
  public:
    /*!
     *  @param system the system that will be measured, this is needed to get
     *  the mean field. It needs to have a member CalculateMeanField that
     *  returns a state_type.
     *  @param one_step_before the phases of the mean field one timestep before the
     *  integration.
     *  @param two_steps_before the  phases of the mean field two timesteps before the
     *  integration.
     *  @param dt the timestep of the integration
     *  @param average_frequency the container, in which the frequency will be
     *  saved.
     */
    AverageFrequencyMeanFieldObserver(system_type& system, 
        const frequency_type& one_step_before, 
        const frequency_type& two_steps_before, double dt, 
        frequency_type& average_frequency)
      : AverageFrequencyObserver<system_type, frequency_type> (system, one_step_before,
          two_steps_before, dt, 0, average_frequency) {}


  protected:
    MeanFieldHelper<mean_field_type> _mean_field_helper;

    frequency_type GetPhases() {
      mean_field_type mean_field = this->_system.CalculateMeanFieldSpherical();
      return _mean_field_helper.GetMeanFieldPhase(mean_field);
    }
};


/*!
 * This observer measures the average frequency of the oscillators and the mean
 * field of a system of phase oscillators. It can be used for general systems as
 * well as networks. This observer is a wrapper around
 * AverageFrequencyPhaseObserver and AverageFrequencyMeanFieldPhaseObserver.
 */
template<typename system_type, typename state_type = std::vector<double>, 
  typename mean_field_type = state_type>
class AverageFrequencyMeanFieldAndPhaseObserver {
  public:
    
    /*!
     *  Detailed constructor, where the user provides the state before the begin
     *  of integration.
     *
     *  @param system the system that will be observed.
     *  @param one_step_before the phases of system state one timestep before the
     *  integration.
     *  @param two_steps_before the phases of the system state two timesteps before the
     *  integration.
     *  @param mean_field_one_before the phases of the mean field one timestep before the
     *  integration.
     *  @param mean_field_two_before the phases of the  mean field two timesteps before the
     *  integration.
     *  @param dt the timestep of the integration
     *  @param dimension the dimension of the ODE.
     *  @param average_frequency_phase the container for the average frequency
     *  of the oscillators.
     *  @param average_frequency_mean_field the container for the average
     *  frequency of the oscillators.
     */
    AverageFrequencyMeanFieldAndPhaseObserver(system_type& system,
        const frequency_type& one_step_before, const frequency_type& two_steps_before,
        const frequency_type& mean_field_one_before,
        const frequency_type& mean_field_two_before,
        double dt, unsigned int dimension, frequency_type& average_frequency_phase,
        frequency_type& average_frequency_mean_field) {
      InitializePointers(system, one_step_before, two_steps_before, 
          mean_field_one_before, mean_field_two_before, dt, dimension,
          average_frequency_phase, average_frequency_mean_field);
    }


    /*!
     *  A simpler constructor where the system is integrated internally for two
     *  timesteps to get the necessary data for the calculation of the average
     *  frequency.
     *
     *  @param system the system that will be integrated.
     *  @param dt the timestep of the integration
     *  @param dimension the dimension of the system.
     *  @param average_frequency_phase the container for the average frequency
     *  of the oscillators.
     *  @param average_frequency_mean_field the container for the average
     *  frequency of the oscillators.
     */
    AverageFrequencyMeanFieldAndPhaseObserver(system_type& system,
        double dt, unsigned int dimension, frequency_type& average_frequency_phase,
        frequency_type& average_frequency_mean_field) {
      frequency_type two_steps_before = GetPhases(system, dimension);
      frequency_type mean_field_two_before = GetPhasesMeanField(system);
      system.Integrate(dt, 1);
      frequency_type one_step_before = GetPhases(system, dimension);
      frequency_type mean_field_one_before = GetPhasesMeanField(system);
      system.Integrate(dt, 1);
      InitializePointers(system, one_step_before, two_steps_before, 
          mean_field_one_before, mean_field_two_before, dt, dimension, 
          average_frequency_phase, average_frequency_mean_field);
    }


    void operator()(const state_type& x, double t) {
      _phase_observer->operator()(x, t);
      _mean_field_observer->operator()(x,t);
    }


  protected:
    std::shared_ptr<AverageFrequencyMeanFieldObserver<system_type, 
      state_type, mean_field_type> > _mean_field_observer;
    std::shared_ptr<AverageFrequencyObserver<system_type, state_type>> 
      _phase_observer;
    MeanFieldHelper<mean_field_type> _mean_field_helper;

    frequency_type GetPhasesMeanField(system_type& system) {
      mean_field_type mean_field = system.CalculateMeanFieldSpherical();
      return _mean_field_helper.GetMeanFieldPhase(mean_field);
    }
    
    frequency_type GetPhases(system_type& system, unsigned int dimension) {
      state_type x = system.GetPositionSpherical();
      frequency_type phases;
      if (dimension == 1) {
        for (size_t i = 0; i < x.size(); ++i) {
          phases.push_back(x[i]);
        }
      }
      else {
        // TODO: check for correct size
        unsigned int N = x.size()/dimension;
        for (unsigned int i = 0; i < N; ++i) {
          for (unsigned int j = i*dimension+1; j < (i+1)*dimension; ++j) {
            phases.push_back(x[j]);
          }
        }
      }
      return phases;
    }

    // This is refactored from the constructors
    void InitializePointers(system_type& system,
        const frequency_type& one_step_before, const frequency_type& two_steps_before,
        const frequency_type& mean_field_one_before,
        const frequency_type& mean_field_two_before,
        double dt, unsigned int dimension, frequency_type& average_frequency_phase,
        frequency_type& average_frequency_mean_field) {
      _mean_field_observer = std::shared_ptr<AverageFrequencyMeanFieldObserver<
        system_type, state_type, mean_field_type> >(
            new AverageFrequencyMeanFieldObserver<
        system_type, state_type, mean_field_type>(system, mean_field_one_before,
            mean_field_two_before, dt, average_frequency_mean_field));
      _phase_observer = std::shared_ptr<AverageFrequencyObserver<
        system_type, state_type> >
        ( new AverageFrequencyObserver<system_type, state_type>(system,
          one_step_before, two_steps_before, dt, dimension, 
          average_frequency_phase));
    }
};
} // nonlinear_systems

#endif
