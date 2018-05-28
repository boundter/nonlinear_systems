#define BOOST_TEST_MODULE MinimalOrderParameterObserver
#include <boost/test/included/unit_test.hpp>
#include <vector>
#include <cmath>
#include <nonlinear_systems/systems/generic_system.hpp>
#include <nonlinear_systems/systems/generic_network.hpp>
#include <nonlinear_systems/odes/harmonic_oscillator_ode.hpp>
#include <nonlinear_systems/observers/minimal_order_parameter_observer.hpp>

typedef std::vector<double> state_type;
typedef std::vector<state_type> network_type;

using namespace nonlinear_systems;


BOOST_AUTO_TEST_CASE(simple_ode) {
  // HarmonicOscillatorODE is just a dummy
  double params = 1.;
  GenericSystem<HarmonicOscillatorODE, state_type> system(2, 1, &params);
  state_type state_1 = {0., M_PI/2.};
  double order_1 = sqrt(2.)/2.;
  state_type state_2 = {1., 1.};
  state_type state_3 = {M_PI, 0};
  state_type minimal_order(1, 0);
  MinimalOrderParameterObserver<GenericSystem<HarmonicOscillatorODE>, state_type>
    observer(system, minimal_order);
  system.SetPosition(state_1);
  observer(state_1, 0.);
  BOOST_REQUIRE_EQUAL(minimal_order.size(), 1);
  BOOST_CHECK_CLOSE(minimal_order[0], order_1, 0.01);
  system.SetPosition(state_2);
  observer(state_2, 0.);
  BOOST_REQUIRE_EQUAL(minimal_order.size(), 1);
  BOOST_CHECK_CLOSE(minimal_order[0], order_1, 0.01);
  system.SetPosition(state_3);
  observer(state_3, 0.);
  BOOST_REQUIRE_EQUAL(minimal_order.size(), 1);
  BOOST_CHECK_SMALL(minimal_order[0], 0.01);
}


BOOST_AUTO_TEST_CASE(network) {
  double params = 1.;
  std::vector<unsigned int> node_sizes = {2., 2.};
  GenericNetwork<HarmonicOscillatorODE> network(node_sizes, 1, &params);
  state_type state_1 = {0., M_PI/2., 0., M_PI/2.};
  state_type order_1 = {sqrt(2.)/2., sqrt(2.)/2.};
  state_type state_2 = {1., 1., M_PI, 0.};
  state_type order_2 = {order_1[0], 0.};
  state_type state_3 = {M_PI, 0., 1, 1};
  state_type order_3 = {0, 0.};
  state_type minimal_order(2, 0);
  MinimalOrderParameterObserver<GenericNetwork<HarmonicOscillatorODE>, state_type,
    network_type> observer(network, minimal_order);
  network.SetPosition(state_1);
  observer(state_1, 0.);
  BOOST_REQUIRE_EQUAL(minimal_order.size(), 2);
  BOOST_CHECK_CLOSE(minimal_order[0], order_1[0], 0.01);
  BOOST_CHECK_CLOSE(minimal_order[1], order_1[1], 0.01);
  network.SetPosition(state_2);
  observer(state_2, 0.);
  BOOST_REQUIRE_EQUAL(minimal_order.size(), 2);
  BOOST_CHECK_CLOSE(minimal_order[0], order_2[0], 0.01);
  BOOST_CHECK_SMALL(minimal_order[1], 0.01);
  network.SetPosition(state_3);
  observer(state_3, 0.);
  BOOST_REQUIRE_EQUAL(minimal_order.size(), 2);
  BOOST_CHECK_SMALL(minimal_order[0], 0.01);
  BOOST_CHECK_SMALL(minimal_order[1], 0.01);
}
