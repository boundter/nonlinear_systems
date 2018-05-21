#ifndef __INSTANTANEOUS_FREQUENCY_OBSERVER__
#define __INSTANTANEOUS_FREQUENCY_OBSERVER__

#include <cmath>
#include <vector>
#include <nonlinear_systems/misc/frequency_helper.hpp>
#include <nonlinear_systems/misc/helper.hpp>

typedef std::vector<double> frequency_type;

namespace nonlinear_systems {
template<typename system_type, typename state_type = std::vector<double>>
/*!
 *  \brief Observes the instantaneous frequency of all phases of the oscillators 
 *  of the system.
 */
class InstantaneousFrequencyObserver {
  public:
    std::vector<frequency_type>& _instantaneous_frequency;
    std::vector<double>& _t;
    system_type& _system;

    /*!
     *  This constructs an observer of the instantaneous frequency using the two
     *  point method for the derivative. To start the method the phases of the
     *  oscillators one and two steps before the integration have to be passed
     *  as inital conditions, as well as the time of the step before.
     *
     *  @param system the observed system.
     *  @param one_step_before phases one timestep before the integration.
     *  @param time time one timestep before the integration.
     *  @param two_steps_before phases two timestep before the integration.
     *  @param dt the timestep.
     *  @param dimension the dimension of the oscillators.
     *  @param instantaneous_frequency the observed frequency.
     *  @param t the observed time.
     */
    InstantaneousFrequencyObserver(system_type& system,
      const frequency_type& one_step_before, double t_before, 
      const frequency_type& two_steps_before, double dt , unsigned int dimension,
      std::vector<frequency_type>& instantaneous_frequency, std::vector<double>& t)
    : _instantaneous_frequency(instantaneous_frequency), _t(t), _system(system){
      _d = dimension;
      _t_before = t_before;
      double modulo = 2*M_PI;
      double limit_step_size = 1.;
      _frequency_helper = std::shared_ptr<FrequencyHelper<frequency_type>>(
          new FrequencyHelper<frequency_type>(one_step_before, two_steps_before, dt, 
            modulo, limit_step_size));
    }


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
        unsigned int dimension, std::vector<frequency_type>& instantaneous_frequency,
        std::vector<double>& t) 
    : _instantaneous_frequency(instantaneous_frequency), _t(t), _system(system){
      SetInitialState(dt, dimension);
    }


    void operator()(const state_type& x, double t) {
      frequency_type phases = this->GetPhases();
      _instantaneous_frequency.push_back(_frequency_helper->operator()(phases));
      _t.push_back(_t_before);
      _t_before = t;
    }

  protected:
    double _t_before;
    unsigned int _d;
    std::shared_ptr<FrequencyHelper<frequency_type>> _frequency_helper;
    
    /*! 
     *  \brief Minimal constructor for use in derived classes.
     *
     *  This constructor does not properly initialize the Observer. It is used
     *  to do the bare minimum for derived classes. This solves the problem of
     *  the GetPhases: If the constructor calls GetPhases it will use the one of
     *  the base class, so the one defined in InstantaneousFrequencyObserver,
     *  but if the derived class calls SetInitialState to initialize the rest of
     *  the class it will use the GetPhases defined in the derived class.
     *
     *  @param system the observed system.
     *  @param instantaneous_frequency the observed frequency.
     *  @param t the observed time.
     */
    InstantaneousFrequencyObserver(system_type& system, 
        std::vector<frequency_type>& instantaneous_frequency, std::vector<double>& t)
    : _instantaneous_frequency(instantaneous_frequency), _t(t), _system(system)
    {}

    
    /*!
     *  \brief Does the rest of the worker of the miniaml constructor.
     *
     *  This method should be called after the minimal constructor to properly
     *  set the initial conditions of the observer. The reason for this is
     *  explained in the minimal constructor.
     *
     *  @param dt the timestep.
     *  @param dimension the dimension of the system, this is only needed, if
     *  the default GetPhases is used.
     */
    void SetInitialState(double dt, unsigned int dimension) {
      _d = dimension;
      frequency_type two_steps_before = this->GetPhases();
      _system.Integrate(dt, 1);
      frequency_type one_step_before = this->GetPhases();
      _t_before = _system.GetTime();
      _system.Integrate(dt, 1);
      double modulo = 2*M_PI;
      double limit_step_size = 1.;
      _frequency_helper = std::shared_ptr<FrequencyHelper<frequency_type>>(
          new FrequencyHelper<frequency_type>(one_step_before, two_steps_before, dt, 
            modulo, limit_step_size));
    }


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


template<typename system_type, typename state_type = std::vector<double>,
  typename mean_field_type = state_type>
class InstantaneousFrequencyMeanFieldObserver: 
  public InstantaneousFrequencyObserver<system_type, std::vector<double>> {
  public:

    InstantaneousFrequencyMeanFieldObserver(system_type& system,
        double dt, std::vector<frequency_type>& instantaneous_frequency, 
        std::vector<double>& t)
    : InstantaneousFrequencyObserver<system_type, frequency_type>(
        system, instantaneous_frequency, t){ 
      this->SetInitialState(dt, 0);
    } 


    InstantaneousFrequencyMeanFieldObserver(system_type& system,
      const frequency_type& one_step_before, double t_before, 
      const frequency_type& two_steps_before, double dt,
      std::vector<frequency_type>& instantaneous_frequency, std::vector<double>& t)
    :InstantaneousFrequencyObserver<system_type, frequency_type>(
        system, one_step_before, t_before, two_steps_before, dt, 0, instantaneous_frequency, t) {}

  protected:
    MeanFieldHelper<mean_field_type> _mean_field_helper;

    frequency_type GetPhases() {
      mean_field_type mean_field = this->_system.CalculateMeanFieldSpherical();
      return _mean_field_helper.GetMeanFieldPhase(mean_field);
    }
};


template<typename system_type, typename state_type = std::vector<double>,
  typename mean_field_type = state_type>
class InstantaneousFrequencyMeanFieldAndPhaseObserver {
  public:

  InstantaneousFrequencyMeanFieldAndPhaseObserver(system_type& system,
      double dt, unsigned int dimension, 
      std::vector<frequency_type>& instantaneous_frequency_phase,
      std::vector<frequency_type>& instantaneous_frequency_mean_field,
      std::vector<double>& t) {
    frequency_type system_two_before = GetPhases(system, dimension);
    frequency_type mean_field_two_before = GetPhasesMeanField(system);
    system.Integrate(dt, 1);
    frequency_type system_one_before = GetPhases(system, dimension);
    frequency_type mean_field_one_before = GetPhasesMeanField(system);
    double t_before = system.GetTime();
    system.Integrate(dt, 1);
    _phase_observer = std::shared_ptr<
      InstantaneousFrequencyObserver<system_type, state_type>>( 
          new InstantaneousFrequencyObserver<system_type, state_type>(
            system, system_one_before, t_before, system_two_before, dt, 
            dimension, instantaneous_frequency_phase, t));
    _mean_field_observer = std::shared_ptr<
      InstantaneousFrequencyMeanFieldObserver<system_type, state_type, mean_field_type>>( 
          new InstantaneousFrequencyMeanFieldObserver<system_type, state_type, mean_field_type>(
            system, mean_field_one_before, t_before, mean_field_two_before, dt, 
            instantaneous_frequency_mean_field, t_dummy));
  }


  void operator()(const state_type& x, double t) {
    _phase_observer->operator()(x, t);
    _mean_field_observer->operator()(x, t);
  }
  
  protected:
    std::shared_ptr<InstantaneousFrequencyObserver<system_type, state_type>>
      _phase_observer;
    std::shared_ptr<InstantaneousFrequencyMeanFieldObserver<system_type, 
      state_type, mean_field_type>> _mean_field_observer;
    MeanFieldHelper<mean_field_type> _mean_field_helper;
    std::vector<double> t_dummy;

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
};
} // nonlinear_systems
#endif
