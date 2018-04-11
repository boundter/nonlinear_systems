#ifndef __RAYLEIGH_MEAN_FIELD_SYSTEM__
#define __RAYLEIGH_MEAN_FIELD_SYSTEM__

#include <functional>
#include <memory>
#include <random>
#include <vector>
#include <nonlinear_systems/odes/rayleigh_mean_field_ode.hpp>
#include <nonlinear_systems/systems/generic_system.hpp>

typedef std::vector<double> state_type;
typedef boost::numeric::odeint::runge_kutta4<state_type> stepper_type;

namespace nonlinear_systems {
/*!
 *  This class describes a system of coupled Rayleigh oscillators. The ODE
 *  for an uncoupled oscillator is
 *
 *  x'' - nonlinearity*(1 - x'^2)*x' + frequency^2*x = 0.
 *
 *  The coupling can be in the X or X'- coordinate, depending on the template
 *  parameter. 
 */
template <typename ode_type = RayleighMeanFieldODEY>
class RayleighMeanFieldSystem
  : public GenericSystem<ode_type, state_type, stepper_type >  {
  public:
    
    /*!
     *  The system gets initailized to a random state in the square [-3, 3] with
     *  normal distributed frequencies with mean 1 and stdev 0.01.
     *  
     *  @param N number of oscillators
     *  @param nonlinearity nonlinearity of the ODE
     *  @param coupling coupling to the mean field
     *  @param seed seed for the rng
     */
    RayleighMeanFieldSystem(unsigned int N, double nonlinearity, 
        double coupling, unsigned long int seed=123456789)
      :GenericSystem<ode_type>(N, 2){
      _rng.seed(seed);
      _nonlinearity = nonlinearity;
      _coupling = coupling;
      this->_N = N;
      double min_x = -3., max_x = 3.;
      SetRandomState(min_x, max_x); 
      _frequency.resize(this->_N);
      double mean_frequency = 1., stdev_frequency = 0.01;
      std::normal_distribution<double> normal(mean_frequency, 
          stdev_frequency);
      std::function<double()> normal_dist = std::bind(normal,
          std::ref(_rng));
      _frequency = SampleDistribution(_frequency.size(), &normal_dist); 
      this->_ode = std::unique_ptr<ode_type>(new ode_type(this->_N, _frequency,
            _nonlinearity, _coupling));
    }

    
    /*!
     *  Returns the frequencies of the oscillators.
     */
    state_type GetFrequency() {
      return _frequency;
    } 

    
    /*!
     *  Sets the frequencies of the oscillators.
     */
    void SetFrequency(state_type new_frequency) {
      if (new_frequency.size() != _frequency.size()) {
        throw std::length_error("New frequency has the wrong length!");
      }
      _frequency = new_frequency;
      this->_ode = std::unique_ptr<ode_type>(new ode_type(this->_N, _frequency, 
            _nonlinearity, _coupling));
    }


    /*!
     *  Sets the positions in the phase space to random points in the square
     *  [min_x, max_x].
     */
    void SetRandomState(double min_x, double max_x) {
        std::uniform_real_distribution<double> uniform(min_x, max_x);
        std::function<double()> uniform_dist = std::bind(uniform,
            std::ref(_rng));
        this->_x = SampleDistribution(this->_x.size(), &uniform_dist);
    }

    
    /*!
     *  Calculates the limit cycle of one oscillator with the nonlinearity of the
     *  system and frequency 1. All oscillators are then placed in distances of
     *  equal time on this limit cycle.
     */
    // TODO: allow change of frequency
    void SetSplayState() {
      RayleighMeanFieldSystem<> one_oscillator(1, _nonlinearity, 0.);
      const double kDT = 1e-2;
      const unsigned int kNTransient = 1e6;
      one_oscillator.Integrate(kDT, kNTransient);
      double T = one_oscillator.CalculateMeanPeriod(5, kDT);
      for (unsigned int i = 0; i < this->_N; ++i) {
        one_oscillator.Integrate(T/static_cast<double>(this->_N+1), 1);
        state_type position = one_oscillator.GetPosition();
        this->_x[2*i] = position[0];
        this->_x[2*i+1] = position[1];
      }
    }
    
    
    /*!
     *  Calculates the period of the system by measuring the times when it
     *  crosses the positive x'-axis.
     *
     *  @param n_average number of crossings to average over
     */
    template <typename observer_type = boost::numeric::odeint::null_observer>
    double CalculateMeanPeriod(unsigned int n_average, double dt, 
        observer_type observer = boost::numeric::odeint::null_observer()) {
      return this->template CalculatePeriod<observer_type>(n_average, dt,
          this->CrossedPositiveYAxis, 1, this->LinearApproximationCrossing, 
          observer);
    }


  protected:
    state_type _frequency;
    std::mt19937_64 _rng;
    double _nonlinearity, _coupling;
    
    
    state_type SampleDistribution(size_t number_samples, 
        std::function<double()>* distribution) {
      state_type samples;
      for (size_t i = 0; i < number_samples; ++i) {
        samples.push_back((*distribution)());
      }
      return samples;
    }

  
    static bool CrossedPositiveYAxis(const state_type& previous_state,
        const state_type& current_state) {
      return ((current_state[1] > 0) and
          (std::signbit(current_state[0]) xor 
           std::signbit(previous_state[0])));   
    }

    
    static double LinearApproximationCrossing(const state_type& previous_state,
        double previous_t, const state_type& current_state, 
        double current_t) {
      return (previous_state[0]*current_t - current_state[0]*previous_t)
        /(previous_state[0] - current_state[0]);
    }
};
} // nonlinear_systems

#endif
