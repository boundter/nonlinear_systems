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


BOOST_AUTO_TEST_CASE(moving_average) {
  state_type x_1 = {0., 1.};
  state_type x_2 = {1., 2.};
  state_type x_3 = x_1;
  state_type analytical_1 = {0., 1.};
  state_type analytical_2 = {1./2., 3./2.};
  state_type analytical_3 = {1./3., 4./3.};
  unsigned int step_number = 1;
  state_type average(2, 0);
  UpdateAverage<std::vector<double> >(x_1, step_number, average);
  BOOST_REQUIRE_EQUAL(step_number, 2);
  BOOST_REQUIRE_EQUAL(average.size(), analytical_1.size());
  BOOST_CHECK_SMALL(average[0], 0.01);
  BOOST_CHECK_CLOSE(average[1], analytical_1[1], 0.01);
  UpdateAverage<std::vector<double> >(x_2, step_number, average);
  BOOST_REQUIRE_EQUAL(step_number, 3);
  BOOST_REQUIRE_EQUAL(average.size(), analytical_2.size());
  BOOST_CHECK_CLOSE(average[0], analytical_2[0], 0.01);
  BOOST_CHECK_CLOSE(average[1], analytical_2[1], 0.01);
  UpdateAverage<std::vector<double> >(x_3, step_number, average);
  BOOST_REQUIRE_EQUAL(step_number, 4);
  BOOST_REQUIRE_EQUAL(average.size(), analytical_3.size());
  BOOST_CHECK_CLOSE(average[0], analytical_3[0], 0.01);
  BOOST_CHECK_CLOSE(average[1], analytical_3[1], 0.01);
}
