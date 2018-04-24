#ifndef __GENERIC_SYSTEM__
#define __GENERIC_SYSTEM__

#include <memory> // unique_ptr
#include <stdexcept> // length_error, 
#include <vector>
#include <boost/numeric/odeint/stepper/runge_kutta4.hpp>
#include <boost/numeric/odeint/integrate/check_adapter.hpp>
#include <boost/numeric/odeint/integrate/integrate_n_steps.hpp>
#include <boost/numeric/odeint/integrate/null_observer.hpp>
#include <nonlinear_systems/misc/coordinates.hpp>

// TODO: Which functions need to be virtual?
// TODO: Use proper error handling
// TODO: Add method to integrate to specific phase
// TODO: Add method to measure period over time
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
      _N = system_size;
      _d = dimension;
      _x.resize(_N*_d);
      _t = 0.;
      _ode = std::unique_ptr<GenericODE>(new GenericODE(parameters)); 
    }


    /*!
     *  \brief Return the position in the state space for all elements.
     */
    state_type GetPosition() {
      return _x;
    }


    /*!
     *  \brief Set the position in the state space.
     */
    void SetPosition(const state_type& new_position) {
      if (new_position.size() != _N*_d) {
        throw std::length_error("Trying to set new position of wrong length!");
      }
      _x = new_position;
    }

    
    /*!
     * \brief Change the number of oscillators.
     */
    void Resize(unsigned int N) {
      _N = N;
      _x.resize(_N*_d);
    }


    /*!
     *  \brief Get the time of the system.
     */
    double GetTime() {
      return _t;
    }


    /*!
     *  \brief Set the parameters for the ODE. This creates a new pointer
     *  to the ODE.
     */
    void SetParameters(void* parameters) {
      _ode = std::unique_ptr<GenericODE>(new GenericODE(parameters));
    }


    /*!
     *  \brief Returns the average position of all elements in the state
     *  space.
     */
    state_type CalculateMeanField() {
      state_type mean_field(_d); 
      for (unsigned int i = 0; i < _N; ++i) {
        for (unsigned int j = 0; j < _d; ++j) {
          mean_field[j] += _x[_d*i+j];
        }
      }
      for (unsigned int i = 0; i < _d; ++i) {
        mean_field[i] /= static_cast<double>(_N);
      }
      return mean_field;
    }


    /*!
     * \breif Calculates the coordinates on a sphere of the same dimension as
     * the phase space. If the deimension is 1, the corrdinates will be wrapped
     * around the unit circle. The first coordinate is the radius and the later
     * ones are the phases. Careful: in 3-d this is not the same as spherical
     * coordinates with polar angle and azimuth!
     */
    state_type CalculateMeanFieldSpherical() {
      // for d = 1 wrap around unit circle
      state_type spherical_mean_field;
      if (_d == 1) {
        double x = 0, y = 0;
        for (size_t i = 0; i < _x.size(); ++i) {
          x += cos(_x[i]);
          y += sin(_x[i]);
        }
        spherical_mean_field.resize(2);
        spherical_mean_field[0] = 1/static_cast<double>(_x.size())*sqrt(x*x + y*y);
        spherical_mean_field[1] = atan2(y, x);
      }
      else {
        spherical_mean_field = CartesianToSpherical(CalculateMeanField());
      }
      return spherical_mean_field;
    }


    /*!
     * \brief Integrate the system whith an observer.
     *
     * The observer is a user-specified struct/class(/function?) that
     * receives the current time and state. If none is specified the
     * null_observer will be used which does nothing.
     *
     * The pointer to member ode needs to be dereferenced here, this will
     * always be to the template parameter GenericODE.
     */ 
    template <typename observer_type = boost::numeric::odeint::null_observer>
      void Integrate(double dt, unsigned int number_steps, 
          observer_type observer = boost::numeric::odeint::null_observer()) {
        _t = boost::numeric::odeint::integrate_n_steps(_stepper, (*_ode),
            _x, _t, dt, number_steps, observer);
      }


    // TODO: Check for NULL-Pointer
    // TODO: remove infinte loop
    /*!
     *  \brief Calculate the period of the average of all elements in the
     *  state space.
     *
     *  This function calculates the period of the mean field. For one
     *  oscillator this is the same as using the position. The period is 
     *  calculated by measuring the time between crossings of the Poincare
     *  manifold, a surface in the state space that is crossed transversally
     *  by the trajectory.
     *
     *  @param CrossedPoincareManifold This is a user-defined funtion. It
     *  recieves the previous state and the current state and should
     *  return True, if the mean field crossed the Poincare manifold which
     *  defines the surface where the period is measured. Otherwise it
     *  should return False.
     *
     *  @param n_crossings Number of crossings of the Poincare manifold
     *  between periods. This is useful for oscillators with higher
     *  winding numbers.
     *
     *  @param ApproximateCrossingPoincareManifold This is a user-defined
     *  function. It receives the previous and current time and state from
     *  this it should approximate the time of crossing. If none is
     *  specified it will approximate the time by taking the middle between
     *  the time of the previous and current step.
     */
    template <typename observer_type = boost::numeric::odeint::null_observer>
    double CalculatePeriod(unsigned int n_average, double dt,
        bool (*CrossedPoincareManifold)(
          const state_type& /*previous_state*/,
          const state_type& /*current_state*/),
        unsigned int n_crossings = 1,
        double (*ApproximateCrossingPoincareManifold)(
          const state_type& /*previous_state*/, double /*previous_t*/,
          const state_type& /*current_state*/, double /*current_t*/)
        = NULL,
        observer_type observer = boost::numeric::odeint::null_observer()) {
      std::vector<double> times_of_crossing;
      state_type previous_state = CalculateMeanField();
      unsigned int n_observed_crossings = 0;
      observer(_x, _t);
      // we need one more time of crossing than periods
      while (times_of_crossing.size() < n_average + 1) {
        Integrate(dt, 1);
        observer(_x, _t);
        state_type current_state = CalculateMeanField();
        if (CrossedPoincareManifold(previous_state, current_state)) {
          n_observed_crossings += 1;
          if (n_observed_crossings == n_crossings) {
            double current_time = GetTime();
            double t_approx;
            if (ApproximateCrossingPoincareManifold == NULL) {
              t_approx = BifurcationZerothOrderCrossingPoincare(
                  previous_state, current_time - dt, current_state,
                  current_time);
            }
            else {
              t_approx = ApproximateCrossingPoincareManifold(
                  previous_state, current_time - dt, current_state, 
                  current_time);
            }
            times_of_crossing.push_back(t_approx);
            n_observed_crossings = 0;
          }
        }
        previous_state = current_state;
      }
      return CalculatePeriodFromCrossings(times_of_crossing);
    }

    
  protected:
    std::unique_ptr<GenericODE> _ode;
    unsigned int _N, _d;
    state_type _x;
    double _t;


    /*!
     *  \brief Initializer for the GenericSystem class. No ODE will be
     *  initialized, it is intended for the use in inherited classes.
     */
    GenericSystem(unsigned int system_size, unsigned int dimension) {
      _N = system_size;
      _d = dimension;
      _x.resize(_N*_d);
      _t = 0.;
    }



  private:
    stepper_type _stepper;


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
