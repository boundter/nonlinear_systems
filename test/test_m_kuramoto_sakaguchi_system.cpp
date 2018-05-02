#define BOOST_TEST_MODULE MKuramotoSakaguchiSystem
#include <boost/test/included/unit_test.hpp>
#include <cmath>
#include <random>
#include <stdexcept>
#include <nonlinear_systems/systems/m_kuramoto_sakaguchi_system.hpp>

typedef std::vector<double> state_type;
typedef std::vector<state_type> network_type;
typedef std::vector<unsigned int> node_size_type;

using namespace nonlinear_systems;

struct F {
  F() {
    node_size = {2, 2};
    frequency = state_type(4);
    state_type coupling_1 = {0.2, 0.4};
    state_type coupling_2 = {0.4, 0.3};
    coupling = {coupling_1, coupling_2};
    state_type phase_shift_1 = {0.5, 0.2};
    state_type phase_shift_2 = {0.4, 0.1};
    phase_shift = {phase_shift_1, phase_shift_2};
  } 

  ~F() {;}

  node_size_type node_size;
  state_type frequency;
  network_type coupling;
  network_type phase_shift;
};


BOOST_FIXTURE_TEST_CASE(explicit_constructor, F) {
  // test with different length
  state_type vector_1(1);
  state_type vector_2(2);
  state_type vector_3(3);
  state_type vector_5(5);

  // Test correct passing of arguments
  BOOST_CHECK_NO_THROW(MKuramotoSakaguchiSystem(frequency, coupling,
        phase_shift, node_size));
  // frequency has wrong length
  BOOST_CHECK_THROW(MKuramotoSakaguchiSystem(vector_3, coupling,
        phase_shift, node_size), std::length_error);
  BOOST_CHECK_THROW(MKuramotoSakaguchiSystem(vector_5, coupling,
        phase_shift, node_size), std::length_error);
  // coupling has wrong length
  BOOST_CHECK_THROW(MKuramotoSakaguchiSystem(frequency, {vector_2, vector_2, vector_2},
        phase_shift, node_size), std::length_error);
  BOOST_CHECK_THROW(MKuramotoSakaguchiSystem(frequency, {vector_2},
        phase_shift, node_size), std::length_error);
  BOOST_CHECK_THROW(MKuramotoSakaguchiSystem(frequency, {vector_3, vector_2},
        phase_shift, node_size), std::length_error);
  BOOST_CHECK_THROW(MKuramotoSakaguchiSystem(frequency, {vector_2, vector_3},
        phase_shift, node_size), std::length_error);
  BOOST_CHECK_THROW(MKuramotoSakaguchiSystem(frequency, {vector_1, vector_2},
        phase_shift, node_size), std::length_error);
  BOOST_CHECK_THROW(MKuramotoSakaguchiSystem(frequency, {vector_2, vector_1},
        phase_shift, node_size), std::length_error);
  // phase shift has wrong length
  BOOST_CHECK_THROW(MKuramotoSakaguchiSystem(frequency, coupling,
        {vector_2, vector_2, vector_2}, node_size), std::length_error);
  BOOST_CHECK_THROW(MKuramotoSakaguchiSystem(frequency, coupling,
        {vector_2}, node_size), std::length_error);
  BOOST_CHECK_THROW(MKuramotoSakaguchiSystem(frequency, coupling,
        {vector_3, vector_2}, node_size), std::length_error);
  BOOST_CHECK_THROW(MKuramotoSakaguchiSystem(frequency, coupling,
        {vector_2, vector_3}, node_size), std::length_error);
  BOOST_CHECK_THROW(MKuramotoSakaguchiSystem(frequency, coupling,
        {vector_1, vector_2}, node_size), std::length_error);
  BOOST_CHECK_THROW(MKuramotoSakaguchiSystem(frequency, coupling,
        {vector_2, vector_1}, node_size), std::length_error);

  // Check correct initialization of phases
  MKuramotoSakaguchiSystem system(frequency, coupling, phase_shift, node_size);
  std::mt19937_64 rng(123456789);
  std::uniform_real_distribution<double> uniform(-M_PI, M_PI);
  state_type phases = system.GetPosition();
  state_type samples;
  for (size_t i = 0; i < frequency.size(); ++i) {
    samples.push_back(uniform(rng));
  }
  BOOST_REQUIRE_EQUAL(phases.size(), samples.size());
  BOOST_CHECK_CLOSE(phases[0], samples[0], 0.01);
  BOOST_CHECK_CLOSE(phases[1], samples[1], 0.01);
  BOOST_CHECK_CLOSE(phases[2], samples[2], 0.01);
  BOOST_CHECK_CLOSE(phases[3], samples[3], 0.01);
}


BOOST_AUTO_TEST_CASE(two_nodes_constructor) {
  double eps = 0.3;
  double freq = 1.;
  node_size_type node_size = {2, 2};
  node_size_type node_size_short(1);
  node_size_type node_size_long(3);

  BOOST_CHECK_NO_THROW(MKuramotoSakaguchiSystem(freq, eps, node_size));
  BOOST_CHECK_THROW(MKuramotoSakaguchiSystem(freq, eps, node_size_short),
      std::length_error);
  BOOST_CHECK_THROW(MKuramotoSakaguchiSystem(freq, eps, node_size_long),
      std::length_error);
  
  // Check correct initialization of phases; also test SetRandomUniformState
  MKuramotoSakaguchiSystem system(freq, eps, node_size);

  std::mt19937_64 rng(123456789);
  std::uniform_real_distribution<double> uniform(-M_PI, M_PI);
  state_type phases = system.GetPosition();
  state_type samples;
  for (size_t i = 0; i < 4; ++i) {
    samples.push_back(uniform(rng));
  }
  BOOST_REQUIRE_EQUAL(phases.size(), samples.size());
  BOOST_CHECK_CLOSE(phases[0], samples[0], 0.01);
  BOOST_CHECK_CLOSE(phases[1], samples[1], 0.01);
  BOOST_CHECK_CLOSE(phases[2], samples[2], 0.01);
  BOOST_CHECK_CLOSE(phases[3], samples[3], 0.01);
}


BOOST_AUTO_TEST_CASE(mean_field) {
  node_size_type node_size = {5, 5};
  double eps = 0.3;
  double freq = 1.;

  MKuramotoSakaguchiSystem system(freq, eps, node_size);

  // Test for correct calculation of mean field
  state_type x = {0., 0., M_PI, M_PI, M_PI/2.,
                  -0.5, -0.3, 1.2, 1.5, 2};
  state_type mean_field_1 = {0.2, M_PI/2};
  state_type mean_field_2 = {0.554315, 0.84};
  system.SetPosition(x);
  network_type measured = system.CalculateMeanField();
  BOOST_REQUIRE_EQUAL(measured.size(), 2);
  BOOST_REQUIRE_EQUAL(measured[0].size(), 2);
  BOOST_REQUIRE_EQUAL(measured[1].size(), 2);
  BOOST_CHECK_CLOSE(measured[0][0], mean_field_1[0], 0.01);
  BOOST_CHECK_CLOSE(measured[0][1], mean_field_1[1], 0.01);
  BOOST_CHECK_CLOSE(measured[1][0], mean_field_2[0], 0.01);
  BOOST_CHECK_CLOSE(measured[1][1], mean_field_2[1], 0.01);

  // Test for generalzied mean field
  state_type generalized_mean_field_1 = {0.6, 0};
  state_type generalized_mean_field_2 = {0.337261, -2.21696958};
  measured = system.CalculateGeneralizedMeanField(2);
  BOOST_REQUIRE_EQUAL(measured.size(), 2);
  BOOST_REQUIRE_EQUAL(measured[0].size(), 2);
  BOOST_REQUIRE_EQUAL(measured[1].size(), 2);
  BOOST_CHECK_CLOSE(measured[0][0], generalized_mean_field_1[0], 0.01);
  BOOST_CHECK_SMALL(measured[0][1], 0.01);
  BOOST_CHECK_CLOSE(measured[1][0], generalized_mean_field_2[0], 0.01);
  BOOST_CHECK_CLOSE(measured[1][1], generalized_mean_field_2[1], 0.01);
}


BOOST_FIXTURE_TEST_CASE(forcing, F) {
  MKuramotoSakaguchiSystem system(frequency, coupling, phase_shift, node_size);
  state_type x = {0., M_PI/2., M_PI/2., M_PI};
  system.SetPosition(x);
  state_type forcing_1 = {0.175814, 2.161817171};
  state_type forcing_2 = {0.200289, 1.715838424};
  network_type measured = system.CalculateForcing();
  BOOST_REQUIRE_EQUAL(measured.size(), 2);
  BOOST_REQUIRE_EQUAL(measured[0].size(), 2);
  BOOST_REQUIRE_EQUAL(measured[1].size(), 2);
  BOOST_CHECK_CLOSE(measured[0][0], forcing_1[0], 0.01);
  BOOST_CHECK_CLOSE(measured[0][1], forcing_1[1], 0.01);
  BOOST_CHECK_CLOSE(measured[1][0], forcing_2[0], 0.01);
  BOOST_CHECK_CLOSE(measured[1][1], forcing_2[1], 0.01);
}


BOOST_AUTO_TEST_CASE(perturbed_clusters) {
  node_size_type node_size_2 = {3, 2};
  state_type freq_2(5);
  state_type zero_vector_2(2, 0.);
  network_type phase_shift_2(2, zero_vector_2);
  network_type coupling_2(2, zero_vector_2);

  MKuramotoSakaguchiSystem system_2(freq_2, coupling_2, phase_shift_2, node_size_2);
  system_2.SetPerturbedClusters(0.5, M_PI);
  state_type x_2 = system_2.GetPosition();
  state_type clusters_2 = {-0.25, 0., 0.25, M_PI-0.25, M_PI+0.25};
  BOOST_REQUIRE_EQUAL(x_2.size(), clusters_2.size());
  BOOST_CHECK_CLOSE(x_2[0], clusters_2[0], 0.01);
  BOOST_CHECK_SMALL(x_2[1], 0.01);
  BOOST_CHECK_CLOSE(x_2[2], clusters_2[2], 0.01);
  BOOST_CHECK_CLOSE(x_2[3], clusters_2[3], 0.01);
  BOOST_CHECK_CLOSE(x_2[4], clusters_2[4], 0.01);
  
  node_size_type node_size_3 = {3, 2, 2};
  state_type freq_3(7);
  state_type zero_vector_3(3, 0.);
  network_type phase_shift_3(3, zero_vector_3);
  network_type coupling_3(3, zero_vector_3);

  MKuramotoSakaguchiSystem system_3(freq_3, coupling_3, phase_shift_3, node_size_3);
  system_3.SetPerturbedClusters(0.5, M_PI/2.);
  state_type x_3 = system_3.GetPosition();
  state_type clusters_3 = {-0.25, 0., 0.25, M_PI/2.-0.25, M_PI/2.+0.25, M_PI-0.25, M_PI+0.25};
  BOOST_REQUIRE_EQUAL(x_3.size(), clusters_3.size());
  BOOST_CHECK_CLOSE(x_3[0], clusters_3[0], 0.01);
  BOOST_CHECK_SMALL(x_3[1], 0.01);
  BOOST_CHECK_CLOSE(x_3[2], clusters_3[2], 0.01);
  BOOST_CHECK_CLOSE(x_3[3], clusters_3[3], 0.01);
  BOOST_CHECK_CLOSE(x_3[4], clusters_3[4], 0.01);
  BOOST_CHECK_CLOSE(x_3[5], clusters_3[5], 0.01);
  BOOST_CHECK_CLOSE(x_3[6], clusters_3[6], 0.01);
}
