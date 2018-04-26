#define BOOST_TEST_MODULE MKuramotoSakaguchiSystem
#include <boost/test/included/unit_test.hpp>
#include <cmath>
#include <random>
#include <stdexcept>
#include <nonlinear_systems/odes/m_kuramoto_sakaguchi_ode.hpp>
#include <nonlinear_systems/systems/m_kuramoto_sakaguchi_system.hpp>

typedef std::vector<double> state_type;
typedef std::vector<state_type> network_type;
typedef std::vector<unsigned int> node_size_type;

using namespace nonlinear_systems;


BOOST_AUTO_TEST_CASE(test_ODE_constructor) {
  node_size_type node_indices = {0, 2, 3, 5};

  state_type frequency_short(4);
  state_type frequency_long(6);
  state_type frequency_right(5);
  
  state_type filler_1(1);
  state_type filler_2(2);
  state_type filler_3(3);

  network_type coupling_short(2);
  network_type coupling_long(4);
  network_type coupling_correct = {filler_3, filler_3, filler_3};
  network_type coupling_wrong = {filler_3, filler_1, filler_1};

  network_type phase_shift_short(2);
  network_type phase_shift_long(4);
  network_type phase_shift_correct = {filler_3, filler_3, filler_3};
  network_type phase_shift_wrong = {filler_2, filler_1, filler_3};

  BOOST_CHECK_NO_THROW(MKuramotoSakaguchiODE(frequency_right, coupling_correct,
        phase_shift_correct, node_indices));
  BOOST_CHECK_THROW(MKuramotoSakaguchiODE(frequency_short, coupling_correct,
        phase_shift_correct, node_indices), std::length_error);
  BOOST_CHECK_THROW(MKuramotoSakaguchiODE(frequency_long, coupling_correct,
        phase_shift_correct, node_indices), std::length_error);
  BOOST_CHECK_THROW(MKuramotoSakaguchiODE(frequency_right, coupling_short,
        phase_shift_correct, node_indices), std::length_error);
  BOOST_CHECK_THROW(MKuramotoSakaguchiODE(frequency_right, coupling_long,
        phase_shift_correct, node_indices), std::length_error);
  BOOST_CHECK_THROW(MKuramotoSakaguchiODE(frequency_right, coupling_wrong,
        phase_shift_correct, node_indices), std::length_error);
  BOOST_CHECK_THROW(MKuramotoSakaguchiODE(frequency_right, coupling_correct,
        phase_shift_long, node_indices), std::length_error);
  BOOST_CHECK_THROW(MKuramotoSakaguchiODE(frequency_right, coupling_correct,
        phase_shift_short, node_indices), std::length_error);
  BOOST_CHECK_THROW(MKuramotoSakaguchiODE(frequency_right, coupling_correct,
        phase_shift_wrong, node_indices), std::length_error);
}


BOOST_AUTO_TEST_CASE(test_ODE_mean_field) {
  state_type frequency(10);
  network_type coupling(2);
  for (size_t i = 0; i < coupling.size(); ++i) {
    coupling[i] = state_type(2);
  }
  network_type phase_shift(2);
  for (size_t i = 0; i < phase_shift.size(); ++i) {
    phase_shift[i] = state_type(2);
  }  
  node_size_type node_indices = {0, 5, 10};

  MKuramotoSakaguchiODE ode(frequency, coupling, phase_shift, node_indices);

  // Test for correct calculation of mean field
  state_type x = {0., 0., M_PI, M_PI, M_PI/2.,
                  -0.5, -0.3, 1.2, 1.5, 2};
  state_type mean_field_1 = {0.2, M_PI/2};
  state_type mean_field_2 = {0.554315, 0.84};
  network_type mean_field = {mean_field_1, mean_field_2};
  network_type measured = ode.CalculateMeanField(x);
  BOOST_CHECK_CLOSE(measured[0][0], mean_field_1[0], 0.01);
  BOOST_CHECK_CLOSE(measured[0][1], mean_field_1[1], 0.01);
  BOOST_CHECK_CLOSE(measured[1][0], mean_field_2[0], 0.01);
  BOOST_CHECK_CLOSE(measured[1][1], mean_field_2[1], 0.01);
}


BOOST_AUTO_TEST_CASE(test_ODE) {
  node_size_type node_indices = {0, 2, 4};
  state_type frequency = {0., 0.5, 1., 2.};
  state_type coupling_1 = {1., -1.};
  state_type coupling_2 = {0.5, 0.3};
  network_type coupling = {coupling_1, coupling_2};
  state_type phase_shift_1 = {0.25, 0.15};
  state_type phase_shift_2 = {-0.15, 0.05};
  network_type phase_shift = {phase_shift_1, phase_shift_2};

  MKuramotoSakaguchiODE ode(frequency, coupling, phase_shift, node_indices);
  state_type state = {0.25, 0.3, M_PI, 0.8};
  state_type deriv = {-0.05, 0.423, 0.916, 1.899};
  state_type ode_result(4);
  ode(state, ode_result, 0.);
  for (size_t i = 0; i < deriv.size(); ++i) {
    BOOST_CHECK_CLOSE(deriv[i], ode_result[i], 1);
  }
}


BOOST_AUTO_TEST_CASE(test_system_default_constructor) {
  node_size_type node_size = {2, 1, 2};
  node_size_type node_indices = {0, 2, 3, 5};

  state_type frequency_short(4);
  state_type frequency_long(6);
  state_type frequency_right(5);
  
  state_type filler_1(1);
  state_type filler_2(2);
  state_type filler_3(3);

  network_type coupling_short(2);
  network_type coupling_long(4);
  network_type coupling_correct = {filler_3, filler_3, filler_3};
  network_type coupling_wrong = {filler_3, filler_1, filler_2};

  network_type phase_shift_short(2);
  network_type phase_shift_long(4);
  network_type phase_shift_correct = {filler_3, filler_3, filler_3};
  network_type phase_shift_wrong = {filler_2, filler_1, filler_3};

  // Test correct passing of arguments
  BOOST_CHECK_NO_THROW(MKuramotoSakaguchiSystem(frequency_right, coupling_correct,
        phase_shift_correct, node_size));
  BOOST_CHECK_THROW(MKuramotoSakaguchiSystem(frequency_short, coupling_correct,
        phase_shift_correct, node_size), std::length_error);
  BOOST_CHECK_THROW(MKuramotoSakaguchiSystem(frequency_long, coupling_correct,
        phase_shift_correct, node_size), std::length_error);
  BOOST_CHECK_THROW(MKuramotoSakaguchiSystem(frequency_right, coupling_short,
        phase_shift_correct, node_size), std::length_error);
  BOOST_CHECK_THROW(MKuramotoSakaguchiSystem(frequency_right, coupling_long,
        phase_shift_correct, node_size), std::length_error);
  BOOST_CHECK_THROW(MKuramotoSakaguchiSystem(frequency_right, coupling_wrong,
        phase_shift_correct, node_size), std::length_error);
  BOOST_CHECK_THROW(MKuramotoSakaguchiSystem(frequency_right, coupling_correct,
        phase_shift_long, node_size), std::length_error);
  BOOST_CHECK_THROW(MKuramotoSakaguchiSystem(frequency_right, coupling_correct,
        phase_shift_short, node_size), std::length_error);
  BOOST_CHECK_THROW(MKuramotoSakaguchiSystem(frequency_right, coupling_correct,
        phase_shift_wrong, node_size), std::length_error);

  // Check correct initialization of phases
  MKuramotoSakaguchiSystem system(frequency_right, coupling_correct, 
      phase_shift_correct, node_size);

  std::mt19937_64 rng(123456789);
  std::uniform_real_distribution<double> uniform(-M_PI, M_PI);
  state_type phases = system.GetPosition();
  for (size_t i = 0; i < 5; ++i) {
    BOOST_CHECK_CLOSE_FRACTION(phases[i], uniform(rng), 0.01);
  }
}


BOOST_AUTO_TEST_CASE(test_system_2_constructor) {
  double eps = 0.3;
  double freq = 1.;
  node_size_type node_size_correct = {5, 4};
  node_size_type node_size_short = {5};
  node_size_type node_size_long = {5, 4, 6};

  BOOST_CHECK_NO_THROW(MKuramotoSakaguchiSystem(freq, eps, node_size_correct));
  BOOST_CHECK_THROW(MKuramotoSakaguchiSystem(freq, eps, node_size_short),
      std::length_error);
  BOOST_CHECK_THROW(MKuramotoSakaguchiSystem(freq, eps, node_size_long),
      std::length_error);
  
  // Check correct initialization of phases; also test SetRandomUniformState
  MKuramotoSakaguchiSystem system(freq, eps, node_size_correct);

  std::mt19937_64 rng(123456789);
  std::uniform_real_distribution<double> uniform(-M_PI, M_PI);
  state_type phases = system.GetPosition();
  for (size_t i = 0; i < 9; ++i) {
    BOOST_CHECK_CLOSE_FRACTION(phases[i], uniform(rng), 0.01);
  }
}

BOOST_AUTO_TEST_CASE(test_system_mean_field) {
  node_size_type node_size = {5, 5};
  double eps = 0.3;
  double freq = 1.;

  MKuramotoSakaguchiSystem system(freq, eps, node_size);

  // Test for correct calculation of mean field
  state_type x = {0., 0., M_PI, M_PI, M_PI/2.,
                  -0.5, -0.3, 1.2, 1.5, 2};
  state_type mean_field_1 = {0.2, M_PI/2};
  state_type mean_field_2 = {0.554315, 0.84};
  system.SetState(x);
  network_type measured = system.CalculateMeanField();
  BOOST_CHECK_CLOSE(measured[0][0], mean_field_1[0], 0.01);
  BOOST_CHECK_CLOSE(measured[0][1], mean_field_1[1], 0.01);
  BOOST_CHECK_CLOSE(measured[1][0], mean_field_2[0], 0.01);
  BOOST_CHECK_CLOSE(measured[1][1], mean_field_2[1], 0.01);

  // Test for generalzied mean field
  state_type generalized_mean_field_1 = {0.6, 0};
  state_type generalized_mean_field_2 = {0.337261, -2.21696958};
  measured = system.CalculateGeneralizedMeanField(2);
  BOOST_CHECK_CLOSE(measured[0][0], generalized_mean_field_1[0], 0.01);
  BOOST_CHECK_SMALL(measured[0][1], 0.01);
  BOOST_CHECK_CLOSE(measured[1][0], generalized_mean_field_2[0], 0.01);
  BOOST_CHECK_CLOSE(measured[1][1], generalized_mean_field_2[1], 0.01);
}


BOOST_AUTO_TEST_CASE(test_system_forcing) {
  node_size_type node_size = {2, 2};
  state_type freq(4);
  state_type phase_shift_1 = {0.5, 0.2};
  state_type phase_shift_2 = {0.4, 0.1};
  network_type phase_shift = {phase_shift_1, phase_shift_2};
  state_type coupling_1 = {0.2, 0.4};
  state_type coupling_2 = {0.4, 0.3};
  network_type coupling = {coupling_1, coupling_2};
  
  MKuramotoSakaguchiSystem system(freq, coupling, phase_shift, node_size);
  state_type x = {0., M_PI/2., M_PI/2., M_PI};
  system.SetState(x);
  state_type forcing_1 = {0.175814, 2.161817171};
  state_type forcing_2 = {0.200289, 1.715838424};
  network_type measured = system.CalculateForcing();
  BOOST_CHECK_CLOSE(measured[0][0], forcing_1[0], 0.01);
  BOOST_CHECK_CLOSE(measured[0][1], forcing_1[1], 0.01);
  BOOST_CHECK_CLOSE(measured[1][0], forcing_2[0], 0.01);
  BOOST_CHECK_CLOSE(measured[1][1], forcing_2[1], 0.01);
}


BOOST_AUTO_TEST_CASE(SetPerturbedClusters) {
  node_size_type node_size_2 = {3, 2};
  state_type freq_2(5);
  state_type zero_vector_2(2, 0.);
  network_type phase_shift_2(2, zero_vector_2);
  network_type coupling_2(2, zero_vector_2);

  MKuramotoSakaguchiSystem system_2(freq_2, coupling_2, phase_shift_2, node_size_2);
  system_2.SetPerturbedClusters(0.5, M_PI);
  state_type x_2 = system_2.GetPosition();
  state_type clusters_2 = {-0.25, 0., 0.25, M_PI-0.25, M_PI+0.25};
  BOOST_TEST(x_2.size() == clusters_2.size());
  BOOST_CHECK_CLOSE_FRACTION(x_2[0], clusters_2[0], 0.01);
  BOOST_CHECK_SMALL(x_2[1], 0.01);
  BOOST_CHECK_CLOSE_FRACTION(x_2[2], clusters_2[2], 0.01);
  BOOST_CHECK_CLOSE_FRACTION(x_2[3], clusters_2[3], 0.01);
  BOOST_CHECK_CLOSE_FRACTION(x_2[4], clusters_2[4], 0.01);
  
  node_size_type node_size_3 = {3, 2, 2};
  state_type freq_3(7);
  state_type zero_vector_3(3, 0.);
  network_type phase_shift_3(3, zero_vector_3);
  network_type coupling_3(3, zero_vector_3);

  MKuramotoSakaguchiSystem system_3(freq_3, coupling_3, phase_shift_3, node_size_3);
  system_3.SetPerturbedClusters(0.5, M_PI/2.);
  state_type x_3 = system_3.GetPosition();
  state_type clusters_3 = {-0.25, 0., 0.25, M_PI/2.-0.25, M_PI/2.+0.25, M_PI-0.25, M_PI+0.25};
  BOOST_TEST(x_3.size() == clusters_3.size());
  BOOST_CHECK_CLOSE_FRACTION(x_3[0], clusters_3[0], 0.01);
  BOOST_CHECK_SMALL(x_3[1], 0.01);
  BOOST_CHECK_CLOSE_FRACTION(x_3[2], clusters_3[2], 0.01);
  BOOST_CHECK_CLOSE_FRACTION(x_3[3], clusters_3[3], 0.01);
  BOOST_CHECK_CLOSE_FRACTION(x_3[4], clusters_3[4], 0.01);
  BOOST_CHECK_CLOSE_FRACTION(x_3[5], clusters_3[5], 0.01);
  BOOST_CHECK_CLOSE_FRACTION(x_3[6], clusters_3[6], 0.01);
  
  
  node_size_type node_size = {3, 2};

  MKuramotoSakaguchiSystem system(1., 0., node_size);
  system.SetPerturbedClusters(0.5, M_PI);
  state_type x = system.GetPosition();
  state_type clusters = {-0.25, 0., 0.25, M_PI-0.25, M_PI+0.25};
  BOOST_TEST(x.size() == clusters.size());
  BOOST_CHECK_CLOSE_FRACTION(x[0], clusters[0], 0.01);
  BOOST_CHECK_SMALL(x[1], 0.01);
  BOOST_CHECK_CLOSE_FRACTION(x[2], clusters[2], 0.01);
  BOOST_CHECK_CLOSE_FRACTION(x[3], clusters[3], 0.01);
  BOOST_CHECK_CLOSE_FRACTION(x[4], clusters[4], 0.01);
}
