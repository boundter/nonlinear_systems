#define BOOST_TEST_MODULE GenericSystem
#include <boost/test/included/unit_test.hpp>
#include <nonlinear_systems/systems/generic_system.hpp>
#include <stdexcept>

typedef std::vector<double> state_type;
typedef boost::numeric::odeint::runge_kutta4<state_type> stepper_type;

using namespace nonlinear_systems;

class HarmonicOscillator {
  public:
    double omega;

    HarmonicOscillator(void* params) {
      double omega = (*(double*)params);
    };

    void operator()(const state_type& x, state_type& dx, const double t) {
      dx[0] = x[1];
      dx[1] = -omega*omega*x[0]; 
    };
};


BOOST_AUTO_TEST_CASE(test_position) {
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
