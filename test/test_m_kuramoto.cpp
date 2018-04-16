#define BOOST_TEST_MODULE MKuramotoSakaguchiSystem
#include <boost/test/included/unit_test.hpp>
#include <cmath>
#include <nonlinear_systems/odes/m_kuramoto_sakaguchi_ode.hpp>

typedef std::vector<double> state_type;
typedef std::vector<state_type> network_type;
typedef std::vector<unsigned int> node_size_type;

using namespace nonlinear_systems;

BOOST_AUTO_TEST_CASE(test_ODE) {
  state_type frequency(10);
  network_type coupling(2);
  for (size_t i = 0; i < coupling.size(); ++i) {
    coupling[i] = state_type(5);
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
