#ifndef __GENERIC_SYSTEM__
#define __GENERIC_SYSTEM__

#include <vector>
#include <boost/numeric/odeint/stepper/runge_kutta4.hpp>
#include <boost/numeric/odeint/integrate/check_adapter.hpp>
#include <boost/numeric/odeint/integrate/integrate_n_steps.hpp>
#include <stdexcept>

// TODO: Check inheritance; is virtual necessary here?
// TODO: Error handling without exception
// TODO: Poincare surface for mean-field
// TODO: Poincare surface multiple crossings
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

          double CalculatePeriod(unsigned int n_average, double dt,
              bool (*CrossedPoincareManifold)(
                const state_type& /*previous_state*/,
                const state_type& /*current_state*/), 
              unsigned int n_crossings = 1);
          double CalculatePeriod(unsigned int n_average, double dt,
              bool (*CrossedPoincareManifold)(
                const state_type& /*previous_state*/,
                const state_type& /*current_state*/),
              double (*ApproximateCrossingPoincareManifold)(
                const state_type& /*previous_state*/, double /*previous_t*/,
                const state_type& /*current_state*/, double /*current_t*/),
              unsigned int n_crossings = 1);
          
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

          static double BifurcationZerothOrderCrossingPoincare(
              const state_type& previous_state, double previous_t,
              const state_type& current_state, double current_t);
          double CalculatePeriodFromCrossings(
              const std::vector<double>& times_of_crossing);
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

template <typename GenericODE, typename state_type, typename stepper_type>
double GenericSystem<GenericODE, state_type, stepper_type>::
CalculatePeriod(unsigned int n_average, double dt,
  bool (*CrossedPoincareManifold)(const state_type& /*previous_state*/,
                                  const state_type& /*current_state*/),
  unsigned int n_crossings) {
  return CalculatePeriod(n_average, dt, CrossedPoincareManifold, 
      BifurcationZerothOrderCrossingPoincare, n_crossings);
}
// TODO: Check for NULL-Pointer
// TODO: remove infinte loop
template <typename GenericODE, typename state_type, typename stepper_type>
double GenericSystem<GenericODE, state_type, stepper_type>::
CalculatePeriod(unsigned int n_average, double dt,
    bool (*CrossedPoincareManifold)(const state_type& /*previous_state*/,
                                    const state_type& /*current_state*/),
    double (*ApproximateCrossingPoincareManifold)(
      const state_type& /*previous_state*/, double /*previous_time*/,
      const state_type& /*current_state*/, double /*current_time*/),
    unsigned int n_crossings){
  std::vector<double> times_of_crossing;
  state_type previous_state = CalculateMeanField();
  unsigned int n_observed_crossings = 0;
  // we need one more time of crossing than periods
  while (times_of_crossing.size() < n_average + 1) {
    Integrate(dt, 1);
    state_type current_state = CalculateMeanField();
    if (CrossedPoincareManifold(previous_state, current_state)) {
      n_observed_crossings += 1;
      if (n_observed_crossings == n_crossings) {
        double current_time = GetTime();
        double t_approx = ApproximateCrossingPoincareManifold(previous_state, 
            current_time - dt, current_state, current_time);
        times_of_crossing.push_back(t_approx);
        n_observed_crossings = 0;
      }
    }
    previous_state = current_state;
  }
  return CalculatePeriodFromCrossings(times_of_crossing);
}


template <typename GenericODE, typename state_type, typename stepper_type>
double GenericSystem<GenericODE, state_type, stepper_type>::
BifurcationZerothOrderCrossingPoincare(const state_type& previous_state,
    double previous_time, const state_type& current_state, double current_time) {
  return (current_time + previous_time)/2.;
}

template <typename GenericODE, typename state_type, typename stepper_type>
double GenericSystem<GenericODE, state_type, stepper_type>::
CalculatePeriodFromCrossings(const std::vector<double>& times_of_crossing) {
  double period = 0.;
  for(size_t i = 1; i < times_of_crossing.size(); ++i) {
    period += times_of_crossing[i] - times_of_crossing[i-1];
  }
  return period/static_cast<double>(times_of_crossing.size()-1);
}
#endif
