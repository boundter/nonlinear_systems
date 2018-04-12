#define BOOST_TEST_MODULE GenericNetwork
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

BOOST_AUTO_TEST_CASE(test_conversion_state_matrix) {
  std::vector<unsigned int> node_sizes = {3, 2};
  GenericNetwork<HarmonicOscillator, int> network(node_sizes, 1, NULL);

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
