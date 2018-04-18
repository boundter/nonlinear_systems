#define BOOST_TEST_MODULE MKuramotoSakaguchiSystem
#include <boost/test/included/unit_test.hpp>
#include <cmath>
#include <stdexcept>
#include <nonlinear_systems/odes/m_kuramoto_sakaguchi_ode.hpp>

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
