#define BOOST_TEST_MODULE GenericNetwork
#include <boost/test/included/unit_test.hpp>
#include <cmath>
#include <vector>
#include <nonlinear_systems/odes/harmonic_oscillator_ode.hpp>
#include <nonlinear_systems/systems/generic_network.hpp>

typedef std::vector<double> state_type;
typedef std::vector<std::vector<double> > network_type;

using namespace nonlinear_systems;

// TODO: Test constructor?

class TwoHarmonicOscillators {
  public:
    double omega;

    TwoHarmonicOscillators(void* params) {
      omega = reinterpret_cast<double*>(params)[0];
    }

    void operator()(const state_type& x, state_type& dx, const double t) {
      dx[0] = x[1];
      dx[1] = -omega*omega*x[0]; 
      dx[2] = x[3];
      dx[3] = -2*omega*2*omega*x[2]; 
    }
};


struct F {
  F() {
    node_sizes = {3, 2};
    params = 1.;
    // Harmonic Oscillator is just a dummy
    network = std::unique_ptr<GenericNetwork<HarmonicOscillatorODE> > (
        new GenericNetwork<HarmonicOscillatorODE>(node_sizes, 2, &params));
  }

  std::vector<unsigned int> node_sizes;
  double params;
  std::unique_ptr<GenericNetwork<HarmonicOscillatorODE> > network;
};  


BOOST_FIXTURE_TEST_CASE(conversion_state, F) {
  // set a new state of the correct length
  state_type x = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  BOOST_CHECK_NO_THROW(network->SetPosition(x));
  
  // state was correctly set
  state_type state = network->GetPosition();
  BOOST_REQUIRE_EQUAL(state.size(), x.size());
  BOOST_CHECK_CLOSE(state[0], x[0], 0.01);
  BOOST_CHECK_CLOSE(state[1], x[1], 0.01);
  BOOST_CHECK_CLOSE(state[2], x[2], 0.01);
  BOOST_CHECK_CLOSE(state[3], x[3], 0.01);
  BOOST_CHECK_CLOSE(state[4], x[4], 0.01);
  BOOST_CHECK_CLOSE(state[5], x[5], 0.01);
  BOOST_CHECK_CLOSE(state[6], x[6], 0.01);
  BOOST_CHECK_CLOSE(state[7], x[7], 0.01);
  BOOST_CHECK_CLOSE(state[8], x[8], 0.01);
  BOOST_CHECK_CLOSE(state[9], x[9], 0.01);
  
  // correct node indices
  std::vector<unsigned int> node_indices = {0, 6, 10};
  std::vector<unsigned int> indices = network->GetNodeIndices();
  BOOST_REQUIRE_EQUAL(indices.size(), node_indices.size());
  BOOST_CHECK_EQUAL(indices[0], node_indices[0]);
  BOOST_CHECK_EQUAL(indices[1], node_indices[1]);
  BOOST_CHECK_EQUAL(indices[2], node_indices[2]);

  // output to nodes works correct
  network_type nodes = network->GetNodes();
  BOOST_REQUIRE_EQUAL(nodes.size(), 2);
  BOOST_REQUIRE_EQUAL(nodes[0].size(), 6);
  BOOST_REQUIRE_EQUAL(nodes[1].size(), 4);
  BOOST_CHECK_CLOSE(nodes[0][0], x[0], 0.01);
  BOOST_CHECK_CLOSE(nodes[0][1], x[1], 0.01);
  BOOST_CHECK_CLOSE(nodes[0][2], x[2], 0.01);
  BOOST_CHECK_CLOSE(nodes[0][3], x[3], 0.01);
  BOOST_CHECK_CLOSE(nodes[0][4], x[4], 0.01);
  BOOST_CHECK_CLOSE(nodes[0][5], x[5], 0.01);
  BOOST_CHECK_CLOSE(nodes[1][0], x[6], 0.01);
  BOOST_CHECK_CLOSE(nodes[1][1], x[7], 0.01);
  BOOST_CHECK_CLOSE(nodes[1][2], x[8], 0.01);
  BOOST_CHECK_CLOSE(nodes[1][3], x[9], 0.01);
}


// TODO: Maybe also test with observer?
BOOST_AUTO_TEST_CASE(integrate) {
  std::vector<unsigned int> node_sizes = {1, 1};
  double params = 1.;
  GenericNetwork<TwoHarmonicOscillators> network(node_sizes, 2, &params);
  
  double dt = 0.01;
  unsigned int steps_increase = 100;
  state_type x = {0., 1., 0., 1.};
  network.SetPosition(x);

  // Integrate increases time correctly
  double t0 = network.GetTime();
  network.Integrate(dt, steps_increase);
  BOOST_CHECK_CLOSE(t0 + static_cast<double>(steps_increase)*dt, 
      network.GetTime(), 0.01);

  // Integrate changes the state correctly
  // the analytic solution for the harmonic oscillator is 
  // x(t)=A*sin(omega*t+phi)
  // for x(0) = 0 and x'(0) = 1 the analytic solution is:
  // x(t) = 1/omega*sin(omega*t)
  double sin_1 = sin(dt*steps_increase);
  double cos_1 = cos(dt*steps_increase);
  double sin_2 = sin(2*dt*steps_increase);
  double cos_2 = cos(2*dt*steps_increase);
  state_type final_state = {sin_1, cos_1, 1/2.*sin_2, cos_2};
  state_type state = network.GetPosition();
  for (size_t i = 0; i < final_state.size(); ++i) {
    BOOST_CHECK_CLOSE(final_state[i], state[i], 0.01);
  }
}


BOOST_FIXTURE_TEST_CASE(mean_field, F) {
  state_type state = {1., 3., 2., 4., 3., 5., 1., 3., 2., 4.};
  state_type mean_field_1 = {2., 4.};
  state_type mean_field_2 = {1.5, 3.5};
  network->SetPosition(state);
  network_type mean_field = network->CalculateMeanField();
  BOOST_REQUIRE_EQUAL(mean_field.size(), 2);
  BOOST_REQUIRE_EQUAL(mean_field[0].size(), 2);
  BOOST_REQUIRE_EQUAL(mean_field[1].size(), 2);
  BOOST_CHECK_CLOSE(mean_field[0][0], mean_field_1[0], 0.01);
  BOOST_CHECK_CLOSE(mean_field[0][1], mean_field_1[1], 0.01);
  BOOST_CHECK_CLOSE(mean_field[1][0], mean_field_2[0], 0.01);
  BOOST_CHECK_CLOSE(mean_field[1][1], mean_field_2[1], 0.01);
}


BOOST_AUTO_TEST_CASE(mean_field_spherical) {
  std::vector<unsigned int> node_sizes = {3, 2};
  double  params = 1.;
  
  // 1-dimensional;  Harmonic Oscillator is just a dummy
  GenericNetwork<HarmonicOscillatorODE> network_1(node_sizes, 1, &params);
  state_type x_1 = {0., M_PI/2., M_PI, 0., M_PI/2.};
  network_1.SetPosition(x_1);
  network_type spherical_1 = network_1.CalculateMeanFieldSpherical();
  BOOST_REQUIRE_EQUAL(spherical_1.size(), 2);
  BOOST_REQUIRE_EQUAL(spherical_1[0].size(), 2);
  BOOST_REQUIRE_EQUAL(spherical_1[1].size(), 2);
  BOOST_CHECK_CLOSE(spherical_1[0][0], 1/3., 0.1);
  BOOST_CHECK_CLOSE(spherical_1[0][1], M_PI/2., 0.1);
  BOOST_CHECK_CLOSE(spherical_1[1][0], 0.7071, 0.1);
  BOOST_CHECK_CLOSE(spherical_1[1][1], M_PI/4., 0.1);

  // 2-dimensional;  Harmonic Oscillator is just a dummy
  node_sizes = {4, 4};
  GenericNetwork<HarmonicOscillatorODE> network_2(node_sizes, 2, &params);
  state_type x_2 = {0., 5., 3., 2., 1., 3., 7., 8., 0., 5., 3., 2., 1., 3., 7., 8.};
  network_2.SetPosition(x_2);
  network_type spherical_2 = network_2.CalculateMeanFieldSpherical();
  BOOST_REQUIRE_EQUAL(spherical_2.size(), 2);
  BOOST_REQUIRE_EQUAL(spherical_2[0].size(), 2);
  BOOST_REQUIRE_EQUAL(spherical_2[1].size(), 2);
  BOOST_CHECK_CLOSE(spherical_2[0][0], 5.27, 0.1);
  BOOST_CHECK_CLOSE(spherical_2[0][1], 1.0222, 0.1);
  BOOST_CHECK_CLOSE(spherical_2[1][0], 5.27, 0.1);
  BOOST_CHECK_CLOSE(spherical_2[1][1], 1.0222, 0.1);
  
  // 3-dimensional; ODE will not be used, HarmonicOscillator as dummy
  node_sizes = {3, 3};
  GenericNetwork<HarmonicOscillatorODE> network_3(node_sizes, 3, &params);
  state_type x_3 = {0., 1., 5., 4., 2., 7., 5., 2., 4., 0., 1., 5., 4., 2., 7., 5., 2., 4.};
  network_3.SetPosition(x_3);
  network_type spherical_3 = network_3.CalculateMeanFieldSpherical();
  BOOST_REQUIRE_EQUAL(spherical_3.size(), 2);
  BOOST_REQUIRE_EQUAL(spherical_3[0].size(), 3);
  BOOST_REQUIRE_EQUAL(spherical_3[1].size(), 3);
  BOOST_CHECK_CLOSE(spherical_3[0][0], 6.342, 0.1);
  BOOST_CHECK_CLOSE(spherical_3[0][1], 1.078, 0.1);
  BOOST_CHECK_CLOSE(spherical_3[0][2], 1.2679, 0.1);
  BOOST_CHECK_CLOSE(spherical_3[1][0], 6.342, 0.1);
  BOOST_CHECK_CLOSE(spherical_3[1][1], 1.078, 0.1);
  BOOST_CHECK_CLOSE(spherical_3[1][2], 1.2679, 0.1);
}
