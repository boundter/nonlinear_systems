#define BOOST_TEST_MODULE GenericNetwork
#include <cmath>
#include <boost/test/included/unit_test.hpp>
#include <nonlinear_systems/systems/generic_network.hpp>

typedef std::vector<double> state_type;
typedef std::vector<std::vector<double> > network_type;

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


BOOST_AUTO_TEST_CASE(test_conversion_state_matrix) {
  std::vector<unsigned int> node_sizes = {3, 2};
  double params = 1.;
  GenericNetwork<HarmonicOscillator, int> network(node_sizes, 1, &params);

  // set a new state
  std::vector<int> x = {1, 2, 3, 4, 5};
  BOOST_CHECK_NO_THROW(network.SetState(x));
  
  // state was correctly set
  std::vector<int> state = network.GetState();
  BOOST_TEST(state == x);
  
  // correct node indices
  std::vector<unsigned int> node_indices = {0, 3, 5};
  BOOST_TEST(node_indices == network.GetNodeIndices());

  // output to nodes works correct
  std::vector<std::vector<int> > nodes = network.GetNodes();
  BOOST_TEST(x[0] == nodes[0][0]);
  BOOST_TEST(x[1] == nodes[0][1]);
  BOOST_TEST(x[2] == nodes[0][2]);
  BOOST_TEST(x[3] == nodes[1][0]);
  BOOST_TEST(x[4] == nodes[1][1]);
}


BOOST_AUTO_TEST_CASE(test_time) {
  std::vector<unsigned int> node_sizes = {3, 2};
  double params = 1.;
  GenericNetwork<HarmonicOscillator, int> network(node_sizes, 1, &params);

  BOOST_CHECK_SMALL(network.GetTime(), 0.01);
}


// TODO: Maybe also test with observer?
BOOST_AUTO_TEST_CASE(test_integrate) {
  std::vector<unsigned int> node_sizes = {1, 1};
  double params = 1.;
  GenericNetwork<TwoHarmonicOscillators> network(node_sizes, 2, &params);
  
  double dt = 0.01;
  unsigned int steps_increase = 100;
  state_type x = {0., 1., 0., 1.};
  network.SetState(x);

  // Integrate increases time correctly
  double t0 = network.GetTime();
  network.Integrate(dt, steps_increase);
  BOOST_CHECK_CLOSE(t0, network.GetTime(), 0.01);

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
  state_type state = network.GetState();
  for (size_t i = 0; i < final_state.size(); ++i) {
    BOOST_CHECK_CLOSE(final_state[i], state[i], 0.01);
  }
}


BOOST_AUTO_TEST_CASE(test_mean_field) {
  std::vector<unsigned int> node_sizes = {3, 2};
  double params = 1;
  GenericNetwork<HarmonicOscillator> network(node_sizes, 2, &params);

  state_type state = {1., 3., 2., 4., 3., 5., 1., 3., 2., 4.};
  state_type mean_field_1 = {2., 4.};
  state_type mean_field_2 = {1.5, 3.5};
  network.SetState(state);
  network_type mean_field = network.CalculateMeanField();
  
  BOOST_CHECK_CLOSE(mean_field[0][0], mean_field_1[0], 0.01);
  BOOST_CHECK_CLOSE(mean_field[0][1], mean_field_1[1], 0.01);
  BOOST_CHECK_CLOSE(mean_field[1][0], mean_field_2[0], 0.01);
  BOOST_CHECK_CLOSE(mean_field[1][1], mean_field_2[1], 0.01);
}
