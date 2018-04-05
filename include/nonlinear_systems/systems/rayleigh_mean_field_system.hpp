#ifndef __RAYLEIGH_MEAN_FIELD_SYSTEM__
#define __RAYLEIGH_MEAN_FIELD_SYSTEM__

#include <vector>
#include <random>
#include <functional>
#include <nonlinear_systems/odes/rayleigh_mean_field_ode.hpp>
#include <nonlinear_systems/systems/generic_system.hpp>

typedef std::vector<double> state_type;

namespace nonlinear_systems {
  template <typename ode_type = RayleighMeanFieldODEY>
  class RayleighMeanFieldSystem
    : public GenericSystem<ode_type, std::vector<double>, 
    boost::numeric::odeint::runge_kutta4<state_type> >  {
      public:
        RayleighMeanFieldSystem(unsigned int N, double nonlinearity, 
            double coupling, unsigned long int seed=123456789)
          :GenericSystem<ode_type>(N, 2){
            rng.seed(seed);

            double min_x = -3., max_x = 3.;
            SetRandomState(min_x, max_x);
            
            frequency.resize(N);
            double mean_frequency = 1., stdev_frequency = 0.01;
            std::normal_distribution<double> normal(mean_frequency, 
                stdev_frequency);
            std::function<double()> normal_dist = std::bind(normal,
                std::ref(rng));
            frequency = SampleDistribution(frequency.size(), &normal_dist);
            
            this->ode = new ode_type(N, frequency, nonlinearity, coupling);
          }


        state_type GetFrequency() {
          return frequency;
        } 


        void SetFrequency(state_type new_frequency) {
          if (new_frequency.size() != frequency.size()) {
            throw std::length_error("New frequency has the wrong length!");
          }
          frequency = new_frequency;
        }


        void SetRandomState(double min_x, double max_x) {
            std::uniform_real_distribution<double> uniform(min_x, max_x);
            std::function<double()> uniform_dist = std::bind(uniform,
                std::ref(rng));
            this->x = SampleDistribution(this->x.size(), &uniform_dist);
        }
        
        
        template <typename observer_type = boost::numeric::odeint::null_observer>
        double CalculateMeanPeriod(unsigned int n_average, double dt, 
            observer_type observer = boost::numeric::odeint::null_observer()) {
          return this->template CalculatePeriod<observer_type>(n_average, dt,
              this->CrossedPositiveYAxis, 1, this->LinearApproximationCrossing, observer);
        }


      protected:
        state_type frequency;
        std::mt19937_64 rng;
        
        using GenericSystem<ode_type>::CalculatePeriod;
        using GenericSystem<ode_type>::SetParameters;

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
