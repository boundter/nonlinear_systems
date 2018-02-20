#ifndef __GENERIC_SYSTEM__
#define __GENERIC_SYSTEM__

#include <vector>
#include <boost/numeric/odeint/stepper/runge_kutta4.hpp>
#include <boost/numeric/odeint/integrate/check_adapter.hpp>
#include <boost/numeric/odeint/integrate/integrate_n_steps.hpp>
#include <boost/numeric/odeint/integrate/null_observer.hpp>
#include <stdexcept>

// TODO: Check inheritance; is virtual necessary here?
// TODO: Error handling without exception
namespace nonlinear_systems {
  template<typename GenericODE,
    typename state_type = std::vector<double>,
    typename stepper_type = boost::numeric::odeint::runge_kutta4<state_type> >
      class GenericSystem {
        public:
          /*!
           *  Initialzer for the GenericSystem class
           *  @param system_size number of elements in the system
           *  @param dimension number of ODEs per element
           *  @param parameters pointer to the parameters, which will be passed
           *  on to the ODE
           */
          GenericSystem(unsigned int system_size, unsigned int dimension,
              void* parameters) {
            N = system_size;
            d = dimension;
            x.resize(N*d);
            t = 0.;
            ode = new GenericODE(parameters); 
          }


          // TODO: Change to protected?
          /*!
           *  Initializer for the GenericSystem class. No ODE will be
           *  initialized, it is intended for the use in inherited classes.
           */
          GenericSystem(unsigned int system_size, unsigned int dimension) {
            N = system_size;
            d = dimension;
            x.resize(N*d);
            t = 0.;
          }


          /*!
           *  \brief Return the position in the state space for all elements.
           */
          state_type GetPosition() {
            return x;
          }


          /*!
           *  \brief Set the position in the state space.
           */
          void SetPosition(state_type& new_position) {
            if (new_position.size() != N*d) {
              throw std::length_error("Trying to set new position of wrong length!");
            }
            x = new_position;
          }


          /*!
           *  \brief Get the time of the system.
           */
          double GetTime() {
            return t;
          }


          /*!
           *  \brief Set the parameters for the ODE. This creates a new pointer
           *  to the ODE.
           */
          void SetParameters(void* parameters) {
            ode = new GenericODE(parameters);
          }


          /*!
           *  \brief Returns the average position of all elements in the state
           *  space.
           */
          state_type CalculateMeanField() {
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


          /*!
           * \brief Integrate the system whith an observer.
           *
           * The observer is a user-specified struct/class(/function?) that
           * receives the current time and state. If none is specified the
           * null_observer will be used which does nothing.
           */ 
          template <typename observer_type = boost::numeric::odeint::null_observer>
            void Integrate(double dt, unsigned int number_steps, 
                observer_type observer = boost::numeric::odeint::null_observer()) {
              t = boost::numeric::odeint::integrate_n_steps(stepper, (*ode),
                  x, t, dt, number_steps, observer);
            }


          // TODO: Add an observer
          // TODO: Check for NULL-Pointer
          // TODO: remove infinte loop
          /*!
           *  \brief Calculate the period of the average of all elements in the
           *  state space.
           *
           *  This function calculates the period of the mean field. For one
           *  oscillator this is the same as using the position.
           *
           *  @param CrossedPoincareManifold This is a user-defined funtion. It
           *  recieves the previous state and the current state and should
           *  return True, if the mean field crossed the Poincare manifold which
           *  defines the surface where the period is measured. Otherwise it
           *  should return False.
           *
           *  @param ApproximateCrossingPoincareManifold This is a user-defined
           *  function. It receives the previous and current time and state from
           *  this it should approximate the time of crossing.
           *
           *  @param n_crossings Number of crossings of the Poincare manifold
           *  between periods. This is useful for oscillators with higher
           *  winding numbers.
           */
          double CalculatePeriod(unsigned int n_average, double dt,
              bool (*CrossedPoincareManifold)(
                const state_type& /*previous_state*/,
                const state_type& /*current_state*/),
              double (*ApproximateCrossingPoincareManifold)(
                const state_type& /*previous_state*/, double /*previous_t*/,
                const state_type& /*current_state*/, double /*current_t*/),
              unsigned int n_crossings = 1) {
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


          /*!
           *  \brief Calculate the period of the average of all elements in the
           *  state space.
           *
           * For a more detailed overview look at the upper function. The time
           * of crossing will be approximated as t_before_crossing + dt/2.
           */
          double CalculatePeriod(unsigned int n_average, double dt,
              bool (*CrossedPoincareManifold)(
                const state_type& /*previous_state*/,
                const state_type& /*current_state*/), 
              unsigned int n_crossings = 1) {
            return CalculatePeriod(n_average, dt, CrossedPoincareManifold, 
                BifurcationZerothOrderCrossingPoincare, n_crossings);
          }


        private:
          GenericODE* ode;
          unsigned int N, d;
          state_type x;
          double t;
          stepper_type stepper;


          /*
           *  Return the time of crossing as t_before_crossing + dt/2.
           */
          static double BifurcationZerothOrderCrossingPoincare(
              const state_type& previous_state, double previous_t,
              const state_type& current_state, double current_t) {
            return (previous_t + current_t)/2.;
          }


          /*
           * Calculates the differences between the elements in a vector and
           * averages this difference vector. This is used to calculate the
           * period from measuring the crossing of a Poincare manifold.
           */
          double CalculatePeriodFromCrossings(
              const std::vector<double>& times_of_crossing) {
            double period = 0.;
            for(size_t i = 1; i < times_of_crossing.size(); ++i) {
              period += times_of_crossing[i] - times_of_crossing[i-1];
            }
            return period/static_cast<double>(times_of_crossing.size()-1);
          }
      };
} // nonlinear_systems
#endif
