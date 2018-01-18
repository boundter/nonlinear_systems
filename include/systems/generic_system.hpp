#ifndef __GENERIC_SYSTEM__
#define __GENERIC_SYSTEM__

#include <vector>
#include <boost/numeric/odeint/stepper/runge_kutta4.hpp>
#include <boost/numeric/odeint/integrate/check_adapter.hpp>
#include <boost/numeric/odeint/integrate/integrate_n_steps.hpp>

namespace nonlinear_systems {
  namespace system {
    template<typename GenericODE,
      typename state_type = std::vector<double>,
      typename stepper_type = boost::numeric::odeint::runge_kutta4<state_type> >
        class GenericSystem {
          public:
            GenericSystem(unsigned int system_size, unsigned int dimension,
                void* parameters) {
              N = system_size;
              d = dimension;
              x.resize(N*d);
              t = 0.;
              ode = new GenericODE(parameters);
            };
            
            state_type GetPosition() {
              return x; 
            };
            
            double GetTime() {
              return t; 
            };

            state_type CalculateMeanField() {
              state_type mean_field(d); 
              for (unsigned int i = 0; i < N; ++i) {
                for (unsigned int j = 0; j < d; ++j) {
                  mean_field[j] += x[d*i+j];
                }
              }
              for (unsigned int i = 0; i < d; ++i) {
                mean_field[i] /= (double)N;
              }
              return mean_field;
            };
            
            void SetPosition(state_type& new_position) {
              // TODO: Add size-check
              x = new_position; 
            };
            
            void Integrate(double dt, unsigned int number_steps) {
              t = boost::numeric::odeint::integrate_n_steps(stepper, (*ode),
                  x, t, dt, number_steps);
            };

            template <typename observer_type>
              void Integrate(double dt, unsigned int number_steps, 
                  observer_type observer) {
              t = boost::numeric::odeint::integrate_n_steps(stepper, (*ode),
                  x, t, dt, number_steps, observer);
              };

            void UpdateParameters(void* parameters) {
              ode = new GenericODE(parameters);
            }

          private:
            GenericODE* ode;
            unsigned int N, d;
            state_type x;
            double t;
            stepper_type stepper;
        }
  } // system
} // nonlinear_systems

#endif
