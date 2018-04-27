#define BOOST_TEST_MODULE MKuramotoSakaguchiODE
#include <boost/test/included/unit_test.hpp>
#include <cmath>
#include <stdexcept>
#include <vector>
#include <nonlinear_systems/odes/m_kuramoto_sakaguchi_ode.hpp>

typedef std::vector<double> state_type;
typedef std::vector<state_type> network_type;
typedef std::vector<unsigned int> node_size_type;

using namespace nonlinear_systems;


struct F {
  F() {
    node_indices = {0, 2, 4};
    frequency = {0., 0.5, 1., 2.};
    state_type coupling_1 = {1., -1.};
    state_type coupling_2 = {0.5, 0.3};
    coupling = {coupling_1, coupling_2};
    state_type phase_shift_1 = {0.25, 0.15};
    state_type phase_shift_2 = {-0.15, 0.05};
    phase_shift = {phase_shift_1, phase_shift_2};
  }

  ~ F() {;}
  
  state_type frequency;
  network_type coupling;
  network_type phase_shift;
  node_size_type node_indices;
};


BOOST_FIXTURE_TEST_CASE(constructor, F) {
  // test with different length
  state_type vector_1(1);
  state_type vector_2(2);
  state_type vector_3(3);
  state_type vector_5(5);
  
  // all correct length
  BOOST_CHECK_NO_THROW(MKuramotoSakaguchiODE(frequency, coupling, phase_shift,
        node_indices));
  // frequency has wrong length
  BOOST_CHECK_THROW(MKuramotoSakaguchiODE(vector_3, coupling, phase_shift,
        node_indices), std::length_error);
  BOOST_CHECK_THROW(MKuramotoSakaguchiODE(vector_5, coupling, phase_shift,
        node_indices), std::length_error);
  // coupling has wrong length
  BOOST_CHECK_THROW(MKuramotoSakaguchiODE(frequency, {vector_2}, 
        phase_shift, node_indices), std::length_error);
  BOOST_CHECK_THROW(MKuramotoSakaguchiODE(frequency, {vector_2, vector_2, vector_2}, 
        phase_shift, node_indices), std::length_error);
  BOOST_CHECK_THROW(MKuramotoSakaguchiODE(frequency, {vector_1, vector_2}, 
        phase_shift, node_indices), std::length_error);
  BOOST_CHECK_THROW(MKuramotoSakaguchiODE(frequency, {vector_2, vector_1}, 
        phase_shift, node_indices), std::length_error);
  BOOST_CHECK_THROW(MKuramotoSakaguchiODE(frequency, {vector_3, vector_2}, 
        phase_shift, node_indices), std::length_error);
  BOOST_CHECK_THROW(MKuramotoSakaguchiODE(frequency, {vector_2, vector_3}, 
        phase_shift, node_indices), std::length_error);
  BOOST_CHECK_THROW(MKuramotoSakaguchiODE(frequency, {vector_1, vector_1}, 
        phase_shift, node_indices), std::length_error);
  BOOST_CHECK_THROW(MKuramotoSakaguchiODE(frequency, {vector_3, vector_3}, 
        phase_shift, node_indices), std::length_error);
  // phase_shift has wrong length
  BOOST_CHECK_THROW(MKuramotoSakaguchiODE(frequency, coupling, {vector_2}, 
        node_indices), std::length_error);
  BOOST_CHECK_THROW(MKuramotoSakaguchiODE(frequency, coupling, 
        {vector_2, vector_2, vector_2}, node_indices), std::length_error);
  BOOST_CHECK_THROW(MKuramotoSakaguchiODE(frequency, coupling, 
        {vector_1, vector_2}, node_indices), std::length_error);
  BOOST_CHECK_THROW(MKuramotoSakaguchiODE(frequency, coupling, 
        {vector_2, vector_1}, node_indices), std::length_error);
  BOOST_CHECK_THROW(MKuramotoSakaguchiODE(frequency, coupling, 
        {vector_3, vector_2}, node_indices), std::length_error);
  BOOST_CHECK_THROW(MKuramotoSakaguchiODE(frequency, coupling, 
        {vector_2, vector_3}, node_indices), std::length_error);
  BOOST_CHECK_THROW(MKuramotoSakaguchiODE(frequency, coupling, 
        {vector_1, vector_1}, node_indices), std::length_error);
  BOOST_CHECK_THROW(MKuramotoSakaguchiODE(frequency, coupling, 
        {vector_3, vector_3}, node_indices), std::length_error);
    
  // test for correct initalization
  MKuramotoSakaguchiODE ode(frequency, coupling, phase_shift, node_indices);
  state_type _frequency = ode._frequency;
  BOOST_REQUIRE_EQUAL(_frequency.size(), frequency.size());
  for (size_t i = 0; i < _frequency.size(); ++i) {
    BOOST_CHECK_CLOSE(_frequency[i], frequency[i], 0.01);
  }
  network_type _coupling = ode._coupling;
  BOOST_REQUIRE_EQUAL(_coupling.size(), coupling.size());
  for (size_t i = 0; i < _coupling.size(); ++i) {
    BOOST_REQUIRE_EQUAL(_coupling[i].size(), coupling[i].size());
    for (size_t j = 0; j < _coupling[i].size(); ++j) {
      BOOST_CHECK_CLOSE(_coupling[i][j], coupling[i][j], 0.01);
    }
  }
  network_type _phase_shift = ode._phase_shift;
  BOOST_REQUIRE_EQUAL(_phase_shift.size(), phase_shift.size());
  for (size_t i = 0; i < _phase_shift.size(); ++i) {
    BOOST_REQUIRE_EQUAL(_phase_shift[i].size(), phase_shift[i].size());
    for (size_t j = 0; j < _phase_shift[i].size(); ++j) {
      BOOST_CHECK_CLOSE(_phase_shift[i][j], phase_shift[i][j], 0.01);
    }
  }
  node_size_type _node_indices = ode._node_indices;
  BOOST_REQUIRE_EQUAL(_node_indices.size(), node_indices.size());
  for (size_t i = 0; i < _node_indices.size(); ++i) {
    BOOST_CHECK_EQUAL(_node_indices[i], node_indices[i]);
  }
}


BOOST_FIXTURE_TEST_CASE(mean_field, F) {
  node_indices = {0, 5, 10};
  frequency = std::vector<double>(10);

  MKuramotoSakaguchiODE ode(frequency, coupling, phase_shift, node_indices);
  state_type x = {0., 0., M_PI, M_PI, M_PI/2.,
                  -0.5, -0.3, 1.2, 1.5, 2.};
  state_type mean_field_0 = {0.2, M_PI/2.};
  state_type mean_field_1 = {0.554315, 0.84};
  network_type measured = ode.CalculateMeanField(x);
  BOOST_REQUIRE_EQUAL(measured.size(), node_indices.size() - 1);
  BOOST_REQUIRE_EQUAL(measured[0].size(), 2);
  BOOST_REQUIRE_EQUAL(measured[1].size(), 2);
  BOOST_CHECK_CLOSE(measured[0][0], mean_field_0[0], 0.01);
  BOOST_CHECK_CLOSE(measured[0][1], mean_field_0[1], 0.01);
  BOOST_CHECK_CLOSE(measured[1][0], mean_field_1[0], 0.01);
  BOOST_CHECK_CLOSE(measured[1][1], mean_field_1[1], 0.01);
}



BOOST_FIXTURE_TEST_CASE(integration, F) {
  MKuramotoSakaguchiODE ode(frequency, coupling, phase_shift, node_indices);
  state_type state = {0.25, 0.3, M_PI, 0.8};
  state_type analytical = {-0.05, 0.423, 0.916, 1.899};
  state_type ode_result(4);
  ode(state, ode_result, 0.);
  BOOST_REQUIRE_EQUAL(ode_result.size(), analytical.size());
  BOOST_CHECK_CLOSE(ode_result[0], analytical[0], 1);
  BOOST_CHECK_CLOSE(ode_result[1], analytical[1], 1);
  BOOST_CHECK_CLOSE(ode_result[2], analytical[2], 1);
  BOOST_CHECK_CLOSE(ode_result[3], analytical[3], 1);
}

