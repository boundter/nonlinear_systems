#define BOOST_TEST_MODULE OrderParameterObserver
#include <boost/test/included/unit_test.hpp>
#include <cmath>
#include <vector>
#include <nonlinear_systems/systems/generic_system.hpp>
#include <nonlinear_systems/systems/generic_network.hpp>
#include <nonlinear_systems/odes/harmonic_oscillator_ode.hpp>
#include <nonlinear_systems/observers/order_parameter_observer.hpp>

typedef std::vector<double> state_type;
typedef std::vector<state_type> network_type;

using namespace nonlinear_systems;


BOOST_AUTO_TEST_CASE(variance_state) {
  // HarmonicOscillatorODE is just a dummy
  double params = 1.;
  GenericSystem<HarmonicOscillatorODE, state_type> system(3, 1, &params);
  state_type state_1 = {0., M_PI, M_PI/2.}; // R = 1/3
  state_type state_2 = {M_PI, M_PI, M_PI}; // R = 1
  state_type state_3 = {0., 0., M_PI}; // R = 1/3
  state_type analytical_mean_1 = {1./3.};
  state_type analytical_mean_2 = {2./3.};
  state_type analytical_mean_3 = {5./9.};
  state_type analytical_var_1 = {0.};
  state_type analytical_var_2 = {2./9.};
  state_type analytical_var_3 = {4./27.};
  state_type average(1, 0);
  state_type variance(1, 0);
  VarianceOrderParameterObserver<GenericSystem<HarmonicOscillatorODE> > observer(system, average, variance);
  // need to set the state because the mean field comes from the system
  system.SetPosition(state_1);
  observer(state_1, 0.);
  BOOST_REQUIRE_EQUAL(average.size(), 1);
  BOOST_REQUIRE_EQUAL(variance.size(), 1);
  BOOST_CHECK_CLOSE(average[0], analytical_mean_1[0], 0.01);
  BOOST_CHECK_SMALL(variance[0], 0.01);
  system.SetPosition(state_2);
  observer(state_2, 0.);
  BOOST_REQUIRE_EQUAL(average.size(), 1);
  BOOST_REQUIRE_EQUAL(variance.size(), 1);
  BOOST_CHECK_CLOSE(average[0], analytical_mean_2[0], 0.01);
  BOOST_CHECK_CLOSE(variance[0], analytical_var_2[0], 0.01);
  system.SetPosition(state_3);
  observer(state_3, 0.);
  BOOST_REQUIRE_EQUAL(average.size(), 1);
  BOOST_REQUIRE_EQUAL(variance.size(), 1);
  BOOST_CHECK_CLOSE(average[0], analytical_mean_3[0], 0.01);
  BOOST_CHECK_CLOSE(variance[0], analytical_var_3[0], 0.01);
}


BOOST_AUTO_TEST_CASE(variance_network) {
  // HarmonicOscillatorODE is just a dummy
  double params = 1.;
  std::vector<unsigned int> node_size = {3, 3};
  GenericNetwork<HarmonicOscillatorODE> network(node_size, 1, &params); 
  state_type network_1 = {0., M_PI, M_PI/2., M_PI, M_PI, 0.}; // R = 1/3, 1/3
  state_type network_2 = {M_PI, M_PI, M_PI, 0., M_PI, M_PI/2.}; // R = 1, 1/3
  state_type network_3 = {0., 0., M_PI, 0., 0., 0.}; // R = 1/3, 1
  state_type analytical_mean_1 = {1./3., 1./3.};
  state_type analytical_mean_2 = {2./3., 1./3.};
  state_type analytical_mean_3 = {5./9., 5./9.};
  state_type analytical_var_1 = {0., 0.};
  state_type analytical_var_2 = {2./9., 0.};
  state_type analytical_var_3 = {4./27., 4./27.};
  state_type average(2, 0);
  state_type variance(2, 0);
  VarianceOrderParameterObserver<GenericNetwork<HarmonicOscillatorODE>, state_type,
   network_type > observer(network, average, variance);
  network.SetPosition(network_1);
  observer(network_1, 0.);
  BOOST_REQUIRE_EQUAL(average.size(), 2);
  BOOST_REQUIRE_EQUAL(variance.size(), 2);
  BOOST_CHECK_CLOSE(average[0], analytical_mean_1[0], 0.01);
  BOOST_CHECK_CLOSE(average[1], analytical_mean_1[1], 0.01);
  BOOST_CHECK_SMALL(variance[0], 0.01);
  BOOST_CHECK_SMALL(variance[1], 0.01);
  network.SetPosition(network_2);
  observer(network_2, 0.);
  BOOST_REQUIRE_EQUAL(average.size(), 2);
  BOOST_REQUIRE_EQUAL(variance.size(), 2);
  BOOST_CHECK_CLOSE(average[0], analytical_mean_2[0], 0.01);
  BOOST_CHECK_CLOSE(average[1], analytical_mean_2[1], 0.01);
  BOOST_CHECK_CLOSE(variance[0], analytical_var_2[0], 0.01);
  BOOST_CHECK_SMALL(variance[1], 0.01);
  network.SetPosition(network_3);
  observer(network_3, 0.);
  BOOST_REQUIRE_EQUAL(average.size(), 2);
  BOOST_REQUIRE_EQUAL(variance.size(), 2);
  BOOST_CHECK_CLOSE(average[0], analytical_mean_3[0], 0.01);
  BOOST_CHECK_CLOSE(average[1], analytical_mean_3[1], 0.01);
  BOOST_CHECK_CLOSE(variance[0], analytical_var_3[0], 0.01);
  BOOST_CHECK_CLOSE(variance[1], analytical_var_3[1], 0.01);
}
