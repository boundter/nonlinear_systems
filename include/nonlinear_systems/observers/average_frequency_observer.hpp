#ifndef __AVERAGE_FREQUENCY_OBSERVER__
#define __AVERAGE_FREQUENCY_OBSERVER__

#include <cmath> // copysign, atan2, fabs
#include <vector>
#include <nonlinear_systems/observers/statistics_observer.hpp>
#include <nonlinear_systems/misc/helper.hpp>
#include <nonlinear_systems/misc/derivative.hpp>

// TODO: Add observer for AverageFrequencyMeanFieldObserver
// TODO: Add observer for AverageFrequencyMeanFieldPhaseAndPhaseObserver
// TODO: Simplify observers by combining observers for phase oscillators and
// limit cycle oscillators
namespace nonlinear_systems {

/*!
 * This observer measures the average frequency of phase oscillators. It
 * calculates this for all oscillators using numerical differentiation.
 */
template <typename state_type = std::vector<double> >
class AverageFrequencyPhaseObserver {
  public:
    
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
        std::vector<double>& average_frequency) {
        _one_step_before = one_step_before;
        _two_steps_before = two_steps_before;
        _dt = dt;
        average_frequency.resize(one_step_before.size());
        _average_observer = std::shared_ptr<AverageObserver<std::vector<double>> >(
            new AverageObserver<std::vector<double> >(average_frequency));
      }


    // TODO: unwrapping is only done for 2pi, add more possibilities
    // the differentiation uses the middlepoint method; the numerical derivative
    // x' at the timepoint n with timestep dt is:
    // x'[n] = (x[n+1] - x[n-1])/(2*dt)
    // the phases are temporarily unwrapped at the crossing from -pi to pi or pi
    // to -pi by subtracting/adding 2pi depending on the state before the
    // crossing
    void operator()(const state_type& x, double t) {
      double modulo = 2*M_PI;
      double biggest_step_size = 1.;
      state_type frequency = TwoPointDerivative(_two_steps_before, x, _dt, 
          modulo, biggest_step_size);
      _average_observer->operator()(frequency, t);
      _two_steps_before = _one_step_before;
      _one_step_before = x;
    }
  
  
  protected:
    state_type _one_step_before;
    state_type _two_steps_before;
    double _dt;
    std::shared_ptr<AverageObserver<std::vector<double>> > _average_observer;
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


// TODO: Add tests for network
// TODO: Add tests for non-phase oscillators
/*!
 * This observer measures the average frequency of the mean field of a system of
 * phase oscillators.  It can be used for general systems as well as networks,
 * but the timestep has to be constant.
 */
template<typename system_type, typename state_type = std::vector<double>,
  typename mean_field_type = state_type>
class AverageFrequencyMeanFieldObserver {
  public:
    system_type& _system;

    /*!
     *  @param system the system that will be measured, this is needed to get
     *  the mean field. It needs to have a member CalculateMeanField that
     *  returns a state_type.
     *  @param one_step_before the mean field one timestep before the
     *  integration.
     *  @param two_steps_before the mean field two timesteps before the
     *  integration.
     *  @param dt the timestep of the integration
     *  @param average_frequency the container, in which the frequency will be
     *  saved.
     */
    AverageFrequencyMeanFieldObserver(system_type& system, 
        const mean_field_type& one_step_before, 
        const mean_field_type& two_steps_before, double dt, 
        state_type& average_frequency)
      :_system(system) {
        std::vector<double> mean_field_phase_before = 
          GetMeanFieldPhase(one_step_before);
        std::vector<double> mean_field_phase_two_before = 
          GetMeanFieldPhase(two_steps_before);
        _phase_observer = 
          new AverageFrequencyPhaseObserver<std::vector<double> >(
              mean_field_phase_before, mean_field_phase_two_before,
              dt, average_frequency);
      }


    void operator()(const state_type& x, double t) {
      mean_field_type mean_field = _system.CalculateMeanFieldSpherical();
      std::vector<double> mean_field_phase = GetMeanFieldPhase(mean_field);
      _phase_observer->operator()(mean_field_phase, t);
    }


  protected:
    AverageFrequencyPhaseObserver<std::vector<double> >* _phase_observer;
    MeanFieldHelper<mean_field_type> _mean_field_helper;

    std::vector<double> GetMeanFieldPhase(const mean_field_type& mean_field) {
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
     *  @param system the system that will be integrated.
     *  @param one_step_before the system state one timestep before the
     *  integration.
     *  @param two_steps_before the system state two timesteps before the
     *  integration.
     *  @param mean_field_one_before the mean field one timestep before the
     *  integration.
     *  @param mean_field_two_before the mean field two timesteps before the
     *  integration.
     *  @param dt the timestep of the integration
     *  @param average_frequency_phase the container for the average frequency
     *  of the oscillators.
     *  @param average_frequency_mean_field the container for the average
     *  frequency of the oscillators.
     */
    AverageFrequencyMeanFieldAndPhaseObserver(system_type& system,
        const state_type& one_step_before, const state_type& two_steps_before,
        const mean_field_type& mean_field_one_before,
        const mean_field_type& mean_field_two_before,
        double dt, state_type& average_frequency_phase,
        state_type& average_frequency_mean_field) {
      InitializePointers(system, one_step_before, two_steps_before, 
          mean_field_one_before, mean_field_two_before, dt, 
          average_frequency_phase, average_frequency_mean_field);
    }


    /*!
     *  A simpler constructor where the system is integrated internally for two
     *  timesteps to get the necessary data for the calculation of the average
     *  frequency.
     *
     *  @param system the system that will be integrated.
     *  @param dt the timestep of the integration
     *  @param average_frequency_phase the container for the average frequency
     *  of the oscillators.
     *  @param average_frequency_mean_field the container for the average
     *  frequency of the oscillators.
     */
    AverageFrequencyMeanFieldAndPhaseObserver(system_type& system,
        double dt, state_type& average_frequency_phase,
        state_type& average_frequency_mean_field) {
      state_type two_steps_before = system.GetPosition();
      mean_field_type mean_field_two_before = system.CalculateMeanFieldSpherical();
      system.Integrate(dt, 1);
      state_type one_step_before = system.GetPosition();
      mean_field_type mean_field_one_before = system.CalculateMeanFieldSpherical();
      system.Integrate(dt, 1);
      InitializePointers(system, one_step_before, two_steps_before, 
          mean_field_one_before, mean_field_two_before, dt, 
          average_frequency_phase, average_frequency_mean_field);
    }


    void operator()(const state_type& x, double t) {
      _phase_observer->operator()(x, t);
      _mean_field_observer->operator()(x,t);
    }


  protected:
    std::shared_ptr<AverageFrequencyMeanFieldObserver<system_type, 
      state_type, mean_field_type> >
      _mean_field_observer;
    std::shared_ptr<AverageFrequencyPhaseObserver<state_type> > _phase_observer;

    // This is refactored from the constructors
    void InitializePointers(system_type& system,
        const state_type& one_step_before, const state_type& two_steps_before,
        const mean_field_type& mean_field_one_before,
        const mean_field_type& mean_field_two_before,
        double dt, state_type& average_frequency_phase,
        state_type& average_frequency_mean_field) {
      _mean_field_observer = std::shared_ptr<AverageFrequencyMeanFieldObserver<
        system_type, state_type, mean_field_type> >(
            new AverageFrequencyMeanFieldObserver<
        system_type, state_type, mean_field_type>(system, mean_field_one_before,
            mean_field_two_before, dt, average_frequency_mean_field));
      _phase_observer = std::shared_ptr<AverageFrequencyPhaseObserver<state_type> >
        ( new AverageFrequencyPhaseObserver<state_type>(
          one_step_before, two_steps_before, dt, average_frequency_phase));
    }
};
} // nonlinear_systems

#endif
