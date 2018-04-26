#define BOOST_TEST_MODULE Statistics
#include <boost/test/included/unit_test.hpp>
#include <functional>
#include <random>
#include <vector>
#include <nonlinear_systems/misc/statistics.hpp>

typedef std::vector<double> state_type;

using namespace nonlinear_systems;

struct F {
  F() {rng.seed(seed);}
  
  ~F(){;}

  unsigned long int seed = 123456789;
  std::mt19937_64 rng;
};


BOOST_FIXTURE_TEST_CASE(SampleUniformRealDistribution, F) {
  std::uniform_real_distribution<double> uniform(-1., 1.);
  std::function<double()> uniform_dist = std::bind(uniform, std::ref(rng));
  state_type direct(10);
  for (size_t i = 0; i < direct.size(); ++i) {
    direct[i] = uniform_dist();
  }
  // reseed to get the same samples
  rng.seed(seed);
  state_type sampled = SampleDistribution<state_type, double>(10, &uniform_dist);
  BOOST_REQUIRE_EQUAL(sampled.size(), direct.size());
  for (size_t i = 0; i < sampled.size(); ++i) {
    BOOST_CHECK_CLOSE(sampled[i], direct[i], 0.01);
  }
}
