#define BOOST_TEST_MODULE GenericSystem
#include <boost/test/included/unit_test.hpp>
#include <cmath>
#include <stdexcept>
#include <vector>
#include <boost/numeric/odeint/stepper/runge_kutta4.hpp>
#include <nonlinear_systems/odes/harmonic_oscillator_ode.hpp>
#include <nonlinear_systems/observers/position_observer.hpp>
#include <nonlinear_systems/systems/generic_system.hpp>

typedef std::vector<double> state_type;
typedef boost::numeric::odeint::runge_kutta4<state_type> stepper_type;


using namespace nonlinear_systems;

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


struct Harmonic {
  Harmonic() {
    params[0] = 2.; 
    N = 1;
    dimension = 2;
    system = std::shared_ptr<GenericSystem<HarmonicOscillatorODE> > (
        new GenericSystem<HarmonicOscillatorODE>(N, dimension, params));
  }

  ~Harmonic() {;}

  double params[1];
  unsigned int N;
  unsigned int dimension;
  std::shared_ptr<GenericSystem<HarmonicOscillatorODE> > system;
};


BOOST_FIXTURE_TEST_CASE(position, Harmonic) {
  // state initializes to zero
  state_type position = system->GetPosition();
  BOOST_REQUIRE_EQUAL(position.size(), dimension);
  BOOST_CHECK_SMALL(position[0], 0.01);
  BOOST_CHECK_SMALL(position[1], 0.01);

  // SetPosition correctly sets the internal variable
  state_type new_state {0.5, 0.1};
  system->SetPosition(new_state);
  position = system->GetPosition();
  BOOST_REQUIRE_EQUAL(position.size(), new_state.size());
  BOOST_CHECK_CLOSE(position[0], new_state[0], 0.01);
  BOOST_CHECK_CLOSE(position[1], new_state[1], 0.01);

  // setting the position to a wrong size returns an error
  state_type too_long = {0.3, 0.1, 6.0};
  state_type too_short = {0.1};
  BOOST_CHECK_THROW(system->SetPosition(too_long), std::length_error);
  BOOST_CHECK_THROW(system->SetPosition(too_short), std::length_error);
}


BOOST_FIXTURE_TEST_CASE(resize_method, Harmonic) {
  system->Resize(2);

  // check if state has correct size
  BOOST_REQUIRE_EQUAL(system->GetPosition().size(), 4);

  // check if new state can be set
  state_type new_state = {1., 2., 3., 4.};
  BOOST_CHECK_NO_THROW(system->SetPosition(new_state));
  state_type position = system->GetPosition();
  BOOST_REQUIRE_EQUAL(position.size(), new_state.size());
  BOOST_CHECK_CLOSE(position[0], new_state[0], 0.01);
  BOOST_CHECK_CLOSE(position[1], new_state[1], 0.01);
  BOOST_CHECK_CLOSE(position[2], new_state[2], 0.01);
  BOOST_CHECK_CLOSE(position[3], new_state[3], 0.01);

  // check false lengths
  state_type too_short(3);
  state_type too_long(5);
  BOOST_CHECK_THROW(system->SetPosition(too_short), std::length_error);
  BOOST_CHECK_THROW(system->SetPosition(too_long), std::length_error);
}


BOOST_FIXTURE_TEST_CASE(integrate, Harmonic){
  // Integrate correctly changes the position
  // Analytic solution is x=A*sin(omega*t + phi)
  // for x(0)=0 and \dot{x}(0)=1 and with omega=2
  // the solution is x = 0.5*sin(omega*t), \dot{x} = cos(omega*t)
  double dt = 0.01;
  unsigned int n = 10000;
  double t = dt*static_cast<double>(n);
  state_type initial_condition {0., 1.};
  system->SetPosition(initial_condition);
  double t_0 = system->GetTime();
  system->Integrate(dt, n);
  state_type analytical = {1./params[0]*sin(params[0]*t), cos(params[0]*t)};
  state_type numerical = system->GetPosition();
  BOOST_REQUIRE_EQUAL(numerical.size(), analytical.size());
  BOOST_CHECK_CLOSE(numerical[0], analytical[0], 0.01);
  BOOST_CHECK_CLOSE(numerical[1], analytical[1], 0.01);
  
  // Integrate increases time correctly
  double t_1 = system->GetTime();
  BOOST_CHECK_CLOSE(t, (t_1 - t_0), 0.000001);
}


BOOST_FIXTURE_TEST_CASE(integrate_with_observer, Harmonic){
  // Integrate correctly changes the position
  // Analytic solution is x=A*sin(omega*t + phi)
  // for x(0)=0 and \dot{x}(0)=1 and with omega=2, A=0.5 and phi=0
  double dt = 0.01;
  unsigned int n = 10;
  state_type initial_condition {0., 1.};
  system->SetPosition(initial_condition);
  std::vector<state_type> position;
  std::vector<double> t;
  system->Integrate(dt, n, PositionObserver<state_type>(position, t));
  BOOST_REQUIRE_EQUAL(position.size(), n+1);
  BOOST_REQUIRE_EQUAL(t.size(), n+1);
  for (size_t i = 0; i < t.size(); ++i) {
    BOOST_CHECK_CLOSE(dt*static_cast<double>(i), t[i], 0.0001);
    state_type analytical = {0.5*sin(params[0]*t[i]), cos(params[0]*t[i])};
    state_type numerical = position[i];
    BOOST_REQUIRE_EQUAL(numerical.size(), analytical.size());
    BOOST_CHECK_CLOSE(numerical[0], analytical[0], 0.01);
    BOOST_CHECK_CLOSE(numerical[1], analytical[1], 0.01);
  }
}


BOOST_FIXTURE_TEST_CASE(set_parameters, Harmonic){
  // Integrate a given time with new parameters and check closeness
  // if everything is like before, but with omega =/= 1, then A = 1/omega
  double dt = 0.01;
  unsigned int n = 10000;
  double t = dt*static_cast<double>(n);
  double new_params = 3.;
  system->SetParameters(&new_params);
  state_type initial_condition {0., 1.};
  system->SetPosition(initial_condition);
  system->Integrate(dt, n);
  state_type analytical = {1./new_params*sin(new_params*t), cos(new_params*t)};
  state_type numerical = system->GetPosition();
  BOOST_REQUIRE_EQUAL(numerical.size(), analytical.size());
  BOOST_CHECK_CLOSE(numerical[0], analytical[0], 0.01);
  BOOST_CHECK_CLOSE(numerical[1], analytical[1], 0.01);
}


BOOST_FIXTURE_TEST_CASE(calculate_mean_field, Harmonic){
  system->Resize(2);
  state_type initial_condition = {0., 5., 1., -2.5};
  state_type analytical =  {0.5, 1.25};
  system->SetPosition(initial_condition);
  state_type numerical = system->CalculateMeanField();
  BOOST_REQUIRE_EQUAL(numerical.size(), analytical.size());
  BOOST_CHECK_CLOSE(numerical[0], analytical[0], 0.01);
  BOOST_CHECK_CLOSE(numerical[1], analytical[1], 0.01);
}


BOOST_AUTO_TEST_CASE(calculate_mean_field_spherical){
  double params = 1.;

  // 1-dimensional; ODE will not be used, HarmonicOscillator as dummy
  GenericSystem<HarmonicOscillatorODE> system_1(5, 1, &params);
  state_type x_1 = {0., 0., 3*M_PI, 3*M_PI, M_PI/2.};
  system_1.SetPosition(x_1);
  state_type spherical_1 = system_1.CalculateMeanFieldSpherical();
  BOOST_REQUIRE_EQUAL(spherical_1.size(), 2);
  BOOST_CHECK_CLOSE(spherical_1[0], 0.2, 0.01);
  BOOST_CHECK_CLOSE(spherical_1[1], M_PI/2., 0.01);
  
  // 2-dimensional; ODE will not be used, HarmonicOscillator as dummy
  GenericSystem<HarmonicOscillatorODE> system_2(4, 2, &params);
  state_type x_2 = {0., 5., 3., 2., 1., 3., 7., 8.};
  system_2.SetPosition(x_2);
  state_type spherical_2 = system_2.CalculateMeanFieldSpherical();
  BOOST_REQUIRE_EQUAL(spherical_2.size(), 2);
  BOOST_CHECK_CLOSE(spherical_2[0], 5.27, 0.1);
  BOOST_CHECK_CLOSE(spherical_2[1], 1.0222, 0.1);
  
  // 3-dimensional; ODE will not be used, HarmonicOscillator as dummy
  GenericSystem<HarmonicOscillatorODE> system_3(3, 3, &params);
  state_type x_3 = {0., 1., 5., 4., 2., 7., 5., 2., 4.};
  system_3.SetPosition(x_3);
  state_type spherical_3 = system_3.CalculateMeanFieldSpherical();
  BOOST_REQUIRE_EQUAL(spherical_3.size(), 3);
  BOOST_CHECK_CLOSE(spherical_3[0], 6.342, 0.1);
  BOOST_CHECK_CLOSE(spherical_3[1], 1.078, 0.1);
  BOOST_CHECK_CLOSE(spherical_3[2], 1.2679, 0.1);
}


BOOST_FIXTURE_TEST_CASE(calculate_period, Harmonic){
  // Integrate the system over 5 periods. For the initial conditions 
  // x(0) = 0 and y(0)=1 with omega = 2 we have T = pi
  double dt = 0.01;
  unsigned int n_average = 5;
  state_type initial_condition {0., 1.};
  system->SetPosition(initial_condition);
  double T = system->CalculatePeriod(n_average, dt, CrossedPositiveYAxis);
  BOOST_CHECK_CLOSE(T, M_PI, 0.1);

  double linear_T = system->CalculatePeriod(n_average, dt, CrossedPositiveYAxis,
      1, LinearApprox);
  BOOST_CHECK_CLOSE(T, M_PI, 0.1);

  // Multiple crossings for a harmonic oscillator should lead to a
  // multiplication of the period
  double double_T = system->CalculatePeriod(n_average, dt, CrossedPositiveYAxis,
      2);
  BOOST_CHECK_CLOSE(double_T, 2*M_PI, 0.1);
}


BOOST_FIXTURE_TEST_CASE(calculate_period_with_observer, Harmonic){
  // Integrate the system over 1 period. For the initial conditions 
  // x(0) = 0 and y(0)=1 with omega = 1 we have T = pi
  double dt = 0.01;
  unsigned int n_average = 1;
  state_type initial_condition {0., 1.};
  system->SetPosition(initial_condition);
  std::vector<state_type> position;
  std::vector<double> t;
  double T = system->CalculatePeriod(n_average, dt, CrossedPositiveYAxis, 1, NULL,
      PositionObserver<state_type>(position, t));
  BOOST_CHECK_CLOSE(T, M_PI, 0.1);
  for (size_t i = 0; i < t.size(); ++i) {
    BOOST_CHECK_CLOSE(dt*static_cast<double>(i), t[i], 0.0001);
    state_type analytical = {0.5*sin(params[0]*t[i]), cos(params[0]*t[i])};
    state_type numerical = position[i];
    BOOST_REQUIRE_EQUAL(numerical.size(), analytical.size());
    BOOST_CHECK_CLOSE(numerical[0], analytical[0], 0.01);
    BOOST_CHECK_CLOSE(numerical[1], analytical[1], 0.01);
  }
}
