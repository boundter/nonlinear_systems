#ifndef __M_KURAMOTO_SAKAGUCHI_SYSTEM__
#define __M_KURAMOTO_SAKAGUCHI_SYSTEM__ 

#include <cmath>
#include <functional>
#include <random>
#include <vector>
#include <nonlinear_systems/odes/m_kuramoto_sakaguchi_ode.hpp>
#include <nonlinear_systems/systems/generic_network.hpp>
#include <nonlinear_systems/misc/statistics.hpp>

typedef std::vector<double> state_type;
typedef std::vector<state_type> network_stype;
typedef std::vector<unsigned int> node_size_type;

namespace nonlinear_systems {
class MKuramotoSakaguchiSystem
  : public GenericNetwork<MKuramotoSakaguchiODE, double> {
  public:

    MKuramotoSakaguchiSystem(const state_type& frequency, 
        const network_type& coupling, const network_type& phase_shift,
        const node_size_type& node_size, unsigned int seed=123456789)
      :GenericNetwork<MKuramotoSakaguchi, double>(node_sizes, 1) {
      _rng.seed(seed);
      SetRandomUniformState();
      _coupling = coupling;
      _phase_shift = phase_shift;
      _node_size = node_size;
      this->_ode = std::unqiue_ptr<MKuramotoSakaguchiODE>(
          new MKuramotoSakaguchiODE(frequency, _coupling, _phase_shift, 
            this->_node_indices));
    }


    MKuramotoSakaguchiSystem(double repulsive_excess, double frequency,
        const node_size_type& node_size, unsigned int seed=123456789)
      : GenericNetwork<MKuramotoSakaguchi, double>(node_size, 1) {
      _rng.seed(seed); 
      SetRandomUniformState();
      state_type frequency_vector(this->_node_indices.back());
      for (size_t i = this->_node_indices[1], i < frequency_vector.size(); ++i) {
        frequency_vector[i] = frequency;   
      }
      state_type coupling_row = {1., -(1.+repulsive_excess)};
      network_type coupling = {coupling_row, coupling_row};
      state_type zero_row = {0., 0.};
      network_type phase_shift = {zero_row, zero_row};
      this->_ode = std::unqiue_ptr<MKuramotoSakaguchiODE>(
          new MKuramotoSakaguchiODE(frequency_vector, coupling, phase_shift, 
            this->_node_indices));
      
    }


    void SetRandomUniformState() {
      std::uniform_real_distribution<double> uniform(-M_PI, M_PI);
      std::function<double()> uniform_dist = std::bind(uniform, std::ref(_rng));
      this->_x = SampleDistribution<state_type, double>(this->_x.size(), 
          &uniform_dist);
    }


    network_type CalculateMeanField() {
      return this->_ode->CalculateMeanField(this->_x);
    }


    network_type CalculateGeneralizedMeanField(int fourier_number) {
      state_type new_state(this->_x.size());
      for (size_t i = 0; i < new_state.size(); ++i) {
        new_state[i] = static_cast<double>(fourier_number)*this->_x[i];
      }
      return this->_ode->CalculateMeanField(new_state);
    }


    network CalculateForcing() {
      network forcing(_node_size.size());
      network_type mean_field = CalculateMeanField();
      for (size_t i = 0; i < forcing.size(); ++i) {
        double real_part = 0., imag_part = 0.;
        for (size_t j = 0; j < _node_size.size(); ++j) {
          double coupling_part = _coupling[i][j]*_node_size[k]
            /this->_node_indices.back()*mean_field[j][0];
          real_part += coupling_part*cos(mean_field[j][1]+_phase_shift[i][j]);
          imag_part += coupling_part*sin(mean_field[j][1]+_phase_shift[i][j]);
        }
        forcing[i][0] = sqrt(real_part*real_part + imag_part*imag_part);
        forcing[i][1] = atan2(imag_part, real_part);
      }
    }


  protected:
    std::mt19937_64 _rng;
    network_type _coupling, _node_size, _phase_shift;
    
}
} // nonlinear_systems
#endif
