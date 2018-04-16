#ifndef __M_KURAMOTO_SAKAGUCHI_ODE__
#define __M_KURAMOTO_SAKAGUCHI_ODE__

#include <cmath> // sqrt, atan2, sin, cos
#include <stdexcept>
#include <vector>

typedef std::vector<double> state_type;
typedef std::vector<state_type> network_type;
typedef std::vector<unsigned int> node_size_type;

namespace nonlinear_systems {
class MKuramotoSakaguchiODE {
  public:
    state_type _frequency;
    network_type _coupling;
    network_type _phase_shift;
    node_size_type _node_indices;


    MKuramotoSakaguchiODE(const state_type& frequency, 
        const network_type& coupling, const network_type& phase_shift,
        const node_size_type& node_indices)
    :_frequency(frequency), _coupling(coupling), _phase_shift(phase_shift),
    _node_indices(node_indices){
      if (_frequency.size() != _node_indices.back()) {
        throw std::length_error("Frequency has the wrong length!");
      }
      if (_coupling.size() != _node_indices.size() - 1) {
        throw std::length_error("Coupling has the wrong first dimension!");
      }
      for (size_t i = 0; i < _coupling.size(); ++i) {
        if (_coupling[i].size() != _node_indices[i+1] - _node_indices[i]) {
          throw std::length_error("Coupling has the wrong second dimension!");
        }
      }
      if (_phase_shift.size() != _node_indices.size() - 1) {
        throw std::length_error("Phase_shift has the wrong dimension!");
      }
      for (size_t i = 0; i < _phase_shift.size(); ++i) {
        if (_phase_shift[i].size() != _node_indices.size() - 1) {
          throw std::length_error("Phase_shift has the wrong second dimension!");
        }
      }
    }


    void operator()(const state_type& x, state_type& dx, const double t) {
      network_type mean_field = CalculateMeanField(x);
      // loop over all nodes (for oscillators)
      for (size_t i = 0; i < _node_indices.size() - 1; ++i) {
        // loop over all oscillators in that node
        for (unsigned int j = _node_indices[i]; j < _node_indices[i+1]; ++j) {
          double sum_coupling = 0.;
          // loop over all nodes for coupling
          for (size_t k = 0; k < mean_field.size(); ++k) {
            double ratio_nodes = static_cast<double>(_node_indices[k+1] - _node_indices[k])
                                 / static_cast<double>(_node_indices.back());
            sum_coupling += _coupling[i][k]*ratio_nodes*mean_field[k][0]
              *sin(mean_field[k][1] - x[j] + _phase_shift[i][k]);
          }
          dx[j] = _frequency[j] + sum_coupling;
        }
      }
    }

    network_type CalculateMeanField(const state_type& x) {
      network_type mean_field;
      for (size_t i = 0; i < _node_indices.size() - 1; ++i) {
        mean_field.push_back(state_type(2));
        unsigned int number_oscillators = _node_indices[i+1] - _node_indices[i];
        double X = 0, Y = 0;
        for (unsigned int j = _node_indices[i]; j < _node_indices[i+1]; ++j) {
          X += cos(x[j]);
          Y += sin(x[j]);
        }
        mean_field[i][0] = 1./static_cast<double>(number_oscillators)
                           *sqrt(X*X + Y*Y);
        mean_field[i][1] = atan2(Y, X);
      }
      return mean_field;
    }
};

}// nonlinear_systems

#endif
