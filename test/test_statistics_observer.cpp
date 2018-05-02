#define BOOST_TEST_MODULE StatisticsObserver
#include <boost/test/included/unit_test.hpp>
#include <vector>
#include <nonlinear_systems/observers/statistics_observer.hpp>

typedef std::vector<double> state_type;

using namespace nonlinear_systems;

struct F {
  F() {
    x_1 = {0., 1.};
    x_2 = {1., 2.};
    x_3 = x_1;
    x_4 = x_2;
    analytical_mean_1 = {0., 1.};
    analytical_var_1 = {0., 0.};
    analytical_mean_2 = {1./2., 3./2.};
    analytical_var_2 = {0.5, 0.5};
    analytical_mean_3 = {1./3., 4./3.};
    analytical_var_3 = {1./3., 1./3.};
    analytical_mean_4 = {2./4., 6./4.};
    analytical_var_4 = {1./3., 1./3.};
  }

  ~F() {;}
 
  state_type x_1, x_2, x_3, x_4;
  state_type analytical_mean_1, analytical_mean_2, analytical_mean_3, analytical_mean_4;
  state_type analytical_var_1, analytical_var_2, analytical_var_3, analytical_var_4;
};


BOOST_FIXTURE_TEST_CASE(average, F) {
  state_type average(2, 0);
  AverageObserver<state_type> observer(average);
  observer(x_1, 0.);
  BOOST_REQUIRE_EQUAL(average.size(), analytical_mean_1.size());
  BOOST_CHECK_SMALL(average[0], 0.01);
  BOOST_CHECK_CLOSE(average[1], analytical_mean_1[1], 0.01);
  observer(x_2, 0.);
  BOOST_REQUIRE_EQUAL(average.size(), analytical_mean_2.size());
  BOOST_CHECK_CLOSE(average[0], analytical_mean_2[0], 0.01);
  BOOST_CHECK_CLOSE(average[1], analytical_mean_2[1], 0.01);
  observer(x_3, 0.);
  BOOST_REQUIRE_EQUAL(average.size(), analytical_mean_3.size());
  BOOST_CHECK_CLOSE(average[0], analytical_mean_3[0], 0.01);
  BOOST_CHECK_CLOSE(average[1], analytical_mean_3[1], 0.01);
  observer(x_4, 0.);
  BOOST_REQUIRE_EQUAL(average.size(), analytical_mean_4.size());
  BOOST_CHECK_CLOSE(average[0], analytical_mean_4[0], 0.01);
  BOOST_CHECK_CLOSE(average[1], analytical_mean_4[1], 0.01);
}


BOOST_FIXTURE_TEST_CASE(variance, F) {
  state_type average(2, 0);
  state_type var(2, 0);
  VarianceObserver<std::vector<double> > observer(average, var);
  observer(x_1, 0.);
  BOOST_REQUIRE_EQUAL(average.size(), analytical_mean_1.size());
  BOOST_REQUIRE_EQUAL(var.size(), analytical_var_1.size());
  BOOST_CHECK_SMALL(average[0], 0.01);
  BOOST_CHECK_CLOSE(average[1], analytical_mean_1[1], 0.01);
  BOOST_CHECK_SMALL(var[0], 0.01);
  BOOST_CHECK_SMALL(var[1], 0.01);
  observer(x_2, 0.);
  BOOST_REQUIRE_EQUAL(average.size(), analytical_mean_2.size());
  BOOST_REQUIRE_EQUAL(var.size(), analytical_var_2.size());
  BOOST_CHECK_CLOSE(average[0], analytical_mean_2[0], 0.01);
  BOOST_CHECK_CLOSE(average[1], analytical_mean_2[1], 0.01);
  BOOST_CHECK_CLOSE(var[0], analytical_var_2[0], 0.01);
  BOOST_CHECK_CLOSE(var[1], analytical_var_2[1], 0.01);
  observer(x_3, 0.);
  BOOST_REQUIRE_EQUAL(average.size(), analytical_mean_3.size());
  BOOST_REQUIRE_EQUAL(var.size(), analytical_var_3.size());
  BOOST_CHECK_CLOSE(average[0], analytical_mean_3[0], 0.01);
  BOOST_CHECK_CLOSE(average[1], analytical_mean_3[1], 0.01);
  BOOST_CHECK_CLOSE(var[0], analytical_var_3[0], 0.01);
  BOOST_CHECK_CLOSE(var[1], analytical_var_3[1], 0.01);
  observer(x_4, 0.);
  BOOST_REQUIRE_EQUAL(average.size(), analytical_mean_4.size());
  BOOST_REQUIRE_EQUAL(var.size(), analytical_var_4.size());
  BOOST_CHECK_CLOSE(average[0], analytical_mean_4[0], 0.01);
  BOOST_CHECK_CLOSE(average[1], analytical_mean_4[1], 0.01);
  BOOST_CHECK_CLOSE(var[0], analytical_var_4[0], 0.01);
  BOOST_CHECK_CLOSE(var[1], analytical_var_4[0], 0.01);
}
