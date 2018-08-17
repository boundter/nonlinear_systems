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


BOOST_AUTO_TEST_CASE(get_name) {
  std::string name = "aa";
  std::string description = "Whatever";
  ArgumentBase base(name, description);
  BOOST_CHECK_EQUAL(base.GetName(), name);
  double dummy_value;
  Argument<double> argument(name, description, dummy_value, 0.0);
  BOOST_CHECK_EQUAL(argument.GetName(), name);
}


BOOST_AUTO_TEST_CASE(get_value_double) {
  double value = 1.5;
  Argument<double> argument("aa", "Whatever", value, 0.);
  std::string value_string = std::to_string(value);
  BOOST_CHECK_EQUAL(argument.GetValueAsString(), value_string);
}


BOOST_AUTO_TEST_CASE(get_value_int) {
  int value = 1;
  Argument<int> argument("aa", "Whatever", value, 0);
  std::string value_string = std::to_string(value);
  BOOST_CHECK_EQUAL(argument.GetValueAsString(), value_string);
}


BOOST_AUTO_TEST_CASE(get_value_string) {
  std::string value = "Test";
  Argument<std::string> argument("aa", "Whatever", value, value);
  std::string value_string = value;
  BOOST_CHECK_EQUAL(argument.GetValueAsString(), value_string);
}


BOOST_AUTO_TEST_CASE(get_value_vector_double) {
  std::vector<double> value = {1.3, 2.5};
  Argument<std::vector<double>> argument("aa", "Whatever", value, value);
  std::string value_string = " 1.3 2.5";
  BOOST_CHECK_EQUAL(argument.GetValueAsString(), value_string);
}


BOOST_AUTO_TEST_CASE(get_value_vector_unsigned_int) {
  std::vector<unsigned int> value = {2, 3};
  Argument<std::vector<unsigned int>> argument("aa", "Whatever", value, value);
  std::string value_string = " 2 3";
  BOOST_CHECK_EQUAL(argument.GetValueAsString(), value_string);
}

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


BOOST_AUTO_TEST_CASE(double_vector_argument) {
  std::vector<double> test_vector;
  int argc = 4;
  char* argv[] = {"test", "--N", "2", "32"};
  std::vector<std::unique_ptr<ArgumentBase>> v;
  v.emplace_back(new Argument<std::vector<double>>("N", "number oscillators",
        test_vector, std::vector<double>()));
  ParseArguments(argc, argv, v);
  BOOST_REQUIRE_EQUAL(test_vector.size(), 2);
  BOOST_CHECK_CLOSE(test_vector[0], 2, 0.01);
  BOOST_CHECK_CLOSE(test_vector[1], 32, 0.01);
}



BOOST_AUTO_TEST_CASE(double_vector_argument_default) {
  std::vector<double> test_vector;
  int argc = 1;
  char* argv[] = {"test"};
  std::vector<std::unique_ptr<ArgumentBase>> v;
  v.emplace_back(new Argument<std::vector<double>>("N", "number oscillators",
        test_vector, {1., 2.}));
  ParseArguments(argc, argv, v);
  BOOST_REQUIRE_EQUAL(test_vector.size(), 2);
  BOOST_CHECK_CLOSE(test_vector[0], 1, 0.01);
  BOOST_CHECK_CLOSE(test_vector[1], 2, 0.01);
}


BOOST_AUTO_TEST_CASE(unsigned_int_vector_argument) {
  std::vector<unsigned int> test_vector;
  int argc = 4;
  char* argv[] = {"test", "--N", "2", "32"};
  std::vector<std::unique_ptr<ArgumentBase>> v;
  v.emplace_back(new Argument<std::vector<unsigned int>>("N", "number oscillators",
        test_vector, std::vector<unsigned int>()));
  ParseArguments(argc, argv, v);
  BOOST_REQUIRE_EQUAL(test_vector.size(), 2);
  BOOST_CHECK_EQUAL(test_vector[0], 2);
  BOOST_CHECK_EQUAL(test_vector[1], 32);
}
