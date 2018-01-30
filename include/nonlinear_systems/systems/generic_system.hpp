#ifndef __GENERIC_SYSTEM__
#define __GENERIC_SYSTEM__

#include <vector>
#include <boost/numeric/odeint/stepper/runge_kutta4.hpp>
#include <boost/numeric/odeint/integrate/check_adapter.hpp>
#include <boost/numeric/odeint/integrate/integrate_n_steps.hpp>
#include <stdexcept>

// TODO: Check inheritance; is virtual necessary here?
// TODO: Error handling without exception
namespace nonlinear_systems {
  template<typename GenericODE,
    typename state_type = std::vector<double>,
    typename stepper_type = boost::numeric::odeint::runge_kutta4<state_type> >
      class GenericSystem {
        public:
          GenericSystem(unsigned int system_size, unsigned int dimension,
              void* parameters);
          state_type GetPosition();
          void SetPosition(state_type& new_position);
          double GetTime();
          void SetParameters(void* parameters);
          state_type CalculateMeanField();
          
          template <typename observer_type>
            void Integrate(double dt, unsigned int number_steps, 
                observer_type observer);
          void Integrate(double dt, unsigned int number_steps);


        private:
          GenericODE* ode;
          unsigned int N, d;
          state_type x;
          double t;
          stepper_type stepper;
      };
} // nonlinear_systems


using namespace nonlinear_systems;


template<typename GenericODE, typename state_type, typename stepper_type>
GenericSystem<GenericODE, state_type, stepper_type>::
GenericSystem(unsigned int system_size, unsigned int dimension, 
    void* parameters) {
  N = system_size;
  d = dimension;
  x.resize(N*d);
  t = 0.;
  ode = new GenericODE(parameters);
}

template<typename GenericODE, typename state_type, typename stepper_type>
state_type GenericSystem<GenericODE, state_type, stepper_type>::
GetPosition() {
  return x;
} 

template<typename GenericODE, typename state_type, typename stepper_type>
double GenericSystem<GenericODE, state_type, stepper_type>::
GetTime() {
  return t;
}

template<typename GenericODE, typename state_type, typename stepper_type>
state_type GenericSystem<GenericODE, state_type, stepper_type>::
CalculateMeanField() {
  state_type mean_field(d); 
  for (unsigned int i = 0; i < N; ++i) {
    for (unsigned int j = 0; j < d; ++j) {
      mean_field[j] += x[d*i+j];

    }
  }
  for (unsigned int i = 0; i < d; ++i) {
    mean_field[i] /= static_cast<double>(N);
  }
  return mean_field;
}

template<typename GenericODE, typename state_type, typename stepper_type>
void GenericSystem<GenericODE, state_type, stepper_type>::
SetPosition(state_type& new_position) {
  if (new_position.size() != N*d) {
    throw std::length_error("Trying to set new position of wrong length!");
  }
  x = new_position; 
}
template<typename GenericODE, typename state_type, typename stepper_type>
void GenericSystem<GenericODE, state_type, stepper_type>::
Integrate(double dt, unsigned int number_steps) {
  t = boost::numeric::odeint::integrate_n_steps(stepper, (*ode),
      x, t, dt, number_steps);
}

template <typename GenericODE, typename state_type, typename stepper_type>
template <typename observer_type>
void GenericSystem<GenericODE, state_type, stepper_type>::
Integrate(double dt, unsigned int number_steps, observer_type observer) {
  t = boost::numeric::odeint::integrate_n_steps(stepper, (*ode),
      x, t, dt, number_steps, observer);
}

template <typename GenericODE, typename state_type, typename stepper_type>
void GenericSystem<GenericODE, state_type, stepper_type>::
SetParameters(void* parameters) {
  ode = new GenericODE(parameters);
}
#endif
