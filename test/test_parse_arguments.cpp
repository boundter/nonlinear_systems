#define BOOST_TEST_MODULE ParseArguments
#include <boost/test/included/unit_test.hpp>
#include <vector>
#include <string>
#include <nonlinear_systems/misc/parse_arguments.hpp>

using namespace nonlinear_systems;


struct F {

  F() {
    v.emplace_back(new Argument<unsigned int>("oscillators", "N", 
          "number oscillators", N, N_0));
    v.emplace_back(new Argument<double>("epsilon", "coupling", eps, eps_0));
    v.emplace_back(new Argument<std::string>("filename", "output", filename, filename_0));
  }

  int argc;
  unsigned int N; 
  unsigned int N_0 = 10;
  double eps;
  double eps_0 = 0.5;
  std::string filename;
  std::string filename_0 = "a.csv";
  std::vector<std::unique_ptr<ArgumentBase>> v;
};

BOOST_FIXTURE_TEST_CASE(default_values, F) {
  int argc = 1;
  char* argv[] = {"test"};
  ParseArguments(argc, argv, v);
  BOOST_CHECK_EQUAL(N, N_0);
  BOOST_CHECK_CLOSE(eps, eps_0, 0.01);
  BOOST_CHECK_EQUAL(filename, filename_0);
}


BOOST_FIXTURE_TEST_CASE(one_argument, F) {
  int argc = 3;
  char* argv[] = {"test", "-N", "50"};
  ParseArguments(argc, argv, v);
  BOOST_CHECK_EQUAL(N, 50);
  BOOST_CHECK_CLOSE(eps, eps_0, 0.01);
  BOOST_CHECK_EQUAL(filename, filename_0);
}


BOOST_FIXTURE_TEST_CASE(multiple_arguments, F) {
  int argc = 7;
  char* argv[] = {"test", "-N", "50", "--filename", "b.csv", "--epsilon", "0.5"};
  ParseArguments(argc, argv, v);
  BOOST_CHECK_EQUAL(N, 50);
  BOOST_CHECK_CLOSE(eps, 0.5, 0.01);
  BOOST_CHECK_EQUAL(filename, "b.csv");
}
