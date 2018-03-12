#ifndef __RAYLEIGH_MEAN_FIELD_SYSTEM__
#define __RAYLEIGH_MEAN_FIELD_SYSTEM__

#include <vector>
#include <random>
#include <functional>
#include <nonlinear_systems/odes/rayleigh_mean_field_ode.hpp>
#include <nonlinear_systems/systems/generic_system.hpp>

typedef std::vector<double> state_type;

namespace nonlinear_systems {
  class RayleighMeanFieldSystem
    : public GenericSystem<RayleighMeanFieldODE,std::vector<double>, 
    boost::numeric::odeint::runge_kutta4<state_type>>  {
      public:
        RayleighMeanFieldSystem(unsigned int N, double nonlinearity, 
            double coupling, const char* coupling_coordinate = "y",
            unsigned long int seed=123456789)
          :GenericSystem(N, 2){
            std::mt19937_64 rng(seed);

            double min_x = -3., max_x = 3.;
            std::uniform_real_distribution<double> uniform(min_x, max_x);
            std::function<double()> uniform_dist = std::bind(uniform,
                std::ref(rng));
            x = SampleDistribution(x.size(), &uniform_dist);
            
            frequency.resize(N);
            double mean_frequency = 1., stdev_frequency = 0.01;
            std::normal_distribution<double> normal(mean_frequency, 
                stdev_frequency);
            std::function<double()> normal_dist = std::bind(normal,
                std::ref(rng));
            frequency = SampleDistribution(frequency.size(), &normal_dist);
            
            _coupling_coordinate = coupling_coordinate;
            if(_coupling_coordinate == "x") {
              ode = new RayleighMeanFieldODEX(N, frequency, nonlinearity, 
                  coupling);
            }
            else if(_coupling_coordinate == "y") {
              ode = new RayleighMeanFieldODEY(N, frequency, nonlinearity, 
                  coupling);
            }
            else {
              throw std::invalid_argument("Did not understand coupling_coordinate. Allowed values are 'x'and 'y'.");
            }
          }


        state_type GetFrequency() {
          return frequency;
        } 

      protected:
        const char* _coupling_coordinate;
        state_type frequency;


        state_type SampleDistribution(size_t number_samples, 
            std::function<double()>* distribution) {
          state_type samples;
          for (size_t i = 0; i < number_samples; ++i) {
            samples.push_back((*distribution)());
          }
          return samples;
        }
    };

} // nonlinear_systems

#endif
