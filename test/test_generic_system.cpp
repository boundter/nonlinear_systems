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

struct PositionAndTimeObserver {
  std::vector<state_type>& pos;
  std::vector<double>& ti;

  PositionAndTimeObserver(std::vector<state_type>& position, 
      std::vector<double>& times)
    :pos(position), ti(times) {}

  void operator() (const state_type& x, double t) {
    pos.push_back(x);
    ti.push_back(t);
  }
};

bool CrossedPositiveYAxis(const state_type& previous_state, 
    const state_type& current_state) {
  return (current_state[1] > 0 and 
      (std::signbit(previous_state[0]) != std::signbit(current_state[0])));
}

double LinearApprox(const state_type& previous_state, double previous_time,
    const state_type& current_state, double current_time) {
  return (previous_state[0]*current_time - current_state[0]*previous_time)
    /(current_state[0] - previous_state[0]);
}


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


BOOST_AUTO_TEST_CASE(test_resize) {
  double params[] = {1.};
  GenericSystem<HarmonicOscillator> system(1, 2, params);

  system.Resize(2);

  // check if state has correct size
  BOOST_TEST(system.GetPosition().size() == 4);

  // check if new state can be set
  state_type new_state {0., 0., 1., 1.};
  BOOST_CHECK_NO_THROW(system.SetPosition(new_state));

  // check false lengths
  state_type too_short {0., 0., 1.};
  state_type too_long {0., 0., 1., 1., 2.};
  BOOST_CHECK_THROW(system.SetPosition(too_short), std::length_error);
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
  BOOST_CHECK_SMALL(dt_time_increase*static_cast<double>(n_time_increase) - 
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
    BOOST_CHECK_SMALL(numeric_solution[i] - analytic_solution[i], 0.0001);
  }
}


BOOST_AUTO_TEST_CASE(test_IntegrateObserver){
  double params = 1.;
  GenericSystem<HarmonicOscillator> system = 
    GenericSystem<HarmonicOscillator>(1, 2, &params);

  // Integrate correctly changes the position
  // Analytic solution is x=A*sin(omega*t + phi)
  // for x(0)=0 and \dot{x}(0)=1 and with omega=1 A=1 and phi=0
  double dt = 0.01;
  unsigned int n = 10;
  state_type initial_condition {0., 1.};
  system.SetPosition(initial_condition);
  std::vector<state_type> position;
  std::vector<double> times;
  system.Integrate(dt, n, PositionAndTimeObserver(position, times));
  for (size_t i = 0; i < times.size(); ++i) {
    BOOST_CHECK_SMALL(dt*static_cast<double>(i) - times[i], 0.0001);
    state_type analytic_solution {sin(times[i]), cos(times[i])};
    state_type numeric_solution = position[i];
    for (size_t j = 0; j < numeric_solution.size(); ++j) {
      BOOST_CHECK_SMALL(numeric_solution[j] - analytic_solution[j], 0.0001);
    }
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
    BOOST_CHECK_SMALL(numeric_solution[i] - analytic_solution[i], 0.0001);
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
    BOOST_CHECK_SMALL(numeric_mean_field[i] - analytic_mean_field[i],
        0.0001);
  }
}


BOOST_AUTO_TEST_CASE(test_CalculatePeriod){
  double params = 1.;
  GenericSystem<HarmonicOscillator> system = 
    GenericSystem<HarmonicOscillator>(1, 2, &params);

  // Integrate the system over 5 periods. For the initial conditions 
  // x(0) = 0 and y(0)=1 with omega = 1 we have T = 2*pi
  double dt = 0.01;
  unsigned int n_average = 5;
  state_type initial_condition {0., 1.};
  system.SetPosition(initial_condition);
  double T = system.CalculatePeriod(n_average, dt, CrossedPositiveYAxis);
  BOOST_CHECK_SMALL(T - 2*M_PI, 0.005);

  double linear_T = system.CalculatePeriod(n_average, dt, CrossedPositiveYAxis,
      1, LinearApprox);
  BOOST_CHECK_SMALL(T - 2*M_PI, 0.005);

  // Multiple crossings for a harmonic oscillator should lead to a
  // multiplication of the period
  double double_T = system.CalculatePeriod(n_average, dt, CrossedPositiveYAxis,
      2);
  BOOST_CHECK_SMALL(double_T - 4*M_PI, 0.005);
}


BOOST_AUTO_TEST_CASE(test_CalculatePeriodObserver){
  double params = 1.;
  GenericSystem<HarmonicOscillator> system = 
    GenericSystem<HarmonicOscillator>(1, 2, &params);

  // Integrate the system over 1 period. For the initial conditions 
  // x(0) = 0 and y(0)=1 with omega = 1 we have T = 2*pi
  double dt = 0.01;
  unsigned int n_average = 1;
  state_type initial_condition {0., 1.};
  system.SetPosition(initial_condition);
  std::vector<state_type> position;
  std::vector<double> times;
  double T = system.CalculatePeriod(n_average, dt, CrossedPositiveYAxis, 1, NULL,
      PositionAndTimeObserver(position, times));
  BOOST_CHECK_SMALL(T - 2*M_PI, 0.005);
  for (size_t i = 0; i < times.size(); ++i) {
    BOOST_CHECK_SMALL(dt*static_cast<double>(i) - times[i], 0.0001);
    state_type analytic_solution {sin(times[i]), cos(times[i])};
    state_type numeric_solution = position[i];
    for (size_t j = 0; j < numeric_solution.size(); ++j) {
      BOOST_CHECK_SMALL(numeric_solution[j] - analytic_solution[j], 0.0001);
    }
  }
}
