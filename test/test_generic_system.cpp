#define BOOST_TEST_MODULE GenericSystem
#include <boost/test/included/unit_test.hpp>
#include <nonlinear_systems/systems/generic_system.hpp>
#include <stdexcept>
#include <cmath>
#include <iostream>

typedef std::vector<double> state_type;
typedef boost::numeric::odeint::runge_kutta4<state_type> stepper_type;

using namespace nonlinear_systems;

class HarmonicOscillator {
  public:
    double omega;

    HarmonicOscillator(void* params) {
      omega = reinterpret_cast<double*>(params)[0];
    }

    void operator()(const state_type& x, state_type& dx, const double t) {
      dx[0] = x[1];
      dx[1] = -omega*omega*x[0]; 
    }
};

class TwoHarmonicOscillators {
  public:
    double* omega;

    TwoHarmonicOscillators(void* params) {
      omega = reinterpret_cast<double*>(params);
    }

    void operator()(const state_type& x, state_type& dx, const double t) {
      dx[0] = x[1];
      dx[1] = -omega[0]*omega[0]*x[0];
      dx[2] = x[3];
      dx[3] = -omega[1]*omega[1]*x[2];
    }
};


BOOST_AUTO_TEST_CASE(test_Position) {
  double params[] = {1.};
  GenericSystem<HarmonicOscillator> system = 
    GenericSystem<HarmonicOscillator>(1, 2, params);

  // state initializes to zero
  state_type zero {0., 0.};
  BOOST_TEST(system.GetPosition() == zero);

  // SetPosition correctly sets the internal variable
  state_type new_state {0.5, 0.1};
  system.SetPosition(new_state);
  BOOST_TEST(system.GetPosition() == new_state);

  // setting the position to a wrong size returns an error
  state_type too_long {0.3, 0.1, 6.0};
  BOOST_CHECK_THROW(system.SetPosition(too_long), std::length_error);
}


BOOST_AUTO_TEST_CASE(test_Integrate){
  double params = 1.;
  GenericSystem<HarmonicOscillator> system = 
    GenericSystem<HarmonicOscillator>(1, 2, &params);

  // Integrate increases time correctly
  double dt_time_increase = 0.01;
  unsigned int n_time_increase = 10000;
  double t_0 = system.GetTime();
  system.Integrate(dt_time_increase, n_time_increase);
  double t_1 = system.GetTime();
  BOOST_CHECK_CLOSE(dt_time_increase*static_cast<double>(n_time_increase), 
      (t_1 - t_0), 0.0001);

  // Integrate correctly changes the position
  // Analytic solution is x=A*sin(omega*t + phi)
  // for x(0)=0 and \dot{x}(0)=1 and with omega=1 A=1 and phi=0
  double dt_position = 0.01;
  unsigned int n_position = 10000;
  double t_position = dt_position*static_cast<double>(n_position);
  state_type initial_condition {0., 1.};
  system.SetPosition(initial_condition);
  system.Integrate(dt_position, n_position);
  state_type analytic_solution {sin(t_position), cos(t_position)};
  state_type numeric_solution = system.GetPosition();
  for (size_t i = 0; i < numeric_solution.size(); ++i) {
    BOOST_CHECK_CLOSE_FRACTION(numeric_solution[i], analytic_solution[i], 0.01);
  }
}

BOOST_AUTO_TEST_CASE(test_Parameters){
  double params = 1.;
  GenericSystem<HarmonicOscillator> system =
    GenericSystem<HarmonicOscillator>(1, 2, &params);

  // Integrate a given time with new parameters and check closeness
  // if everything is like before, but with omega =/= 1, then A = 1/omega
  double dt = 0.01;
  unsigned int n = 10000;
  double t = dt*static_cast<double>(n);
  state_type initial_condition {0., 1.};
  system.SetPosition(initial_condition);
  double new_params = 3.;
  system.SetParameters(&new_params);
  system.Integrate(dt, n);
  state_type analytic_solution {1./new_params*sin(new_params*t), cos(new_params*t)};
  state_type numeric_solution = system.GetPosition();
  for (size_t i = 0; i < numeric_solution.size(); ++i) {
    BOOST_CHECK_CLOSE_FRACTION(numeric_solution[i], analytic_solution[i], 0.01);
  }
}

BOOST_AUTO_TEST_CASE(test_MeanField){
 double params[2] = {1., 3.};
  GenericSystem<TwoHarmonicOscillators> system = 
   GenericSystem<TwoHarmonicOscillators>(2, 2, &params); 

  double dt = 0.01;
  unsigned int n = 10000;
  double t = dt*static_cast<double>(n);
  state_type initial_condition {0., 1., 0., 1.};
  system.SetPosition(initial_condition);
  system.Integrate(dt, n);
  state_type numeric_mean_field = system.CalculateMeanField();
  state_type analytic_mean_field {1/2.*(sin(t) + 1./params[1]*sin(params[1]*t)),
                                  1/2.*(cos(t) + cos(params[1]*t))};
  for (size_t i = 0; i < numeric_mean_field.size(); ++i) {
    BOOST_CHECK_CLOSE_FRACTION(numeric_mean_field[i], analytic_mean_field[i],
        0.01);
  }
}
