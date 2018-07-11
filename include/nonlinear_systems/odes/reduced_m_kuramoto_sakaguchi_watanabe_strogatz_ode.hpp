#ifndef __REDUCED_M_KURAMOTO_SAKAGUCHI_WATANABE_STROGATZ_ODE__
#define __REDUCED_M_KURAMOTO_SAKAGUCHI_WATANABE_STROGATZ_ODE__

#include <cmath> // sin,cos
#include <vector>

typedef std::vector<double> state_type;
typedef std::vector<state_type> network_type;
typedef std::vector<unsigned int> node_size_type;

namespace nonlinear_systems {

class ReducedMKuramotoSakaguchiWatanabeStrogatzODE {
  public:
    network_type _constants;
    state_type _coupling, _phase_shift, _omega;
    unsigned int _N;
    node_size_type _node_indices;


    ReducedMKuramotoSakaguchiWatanabeStrogatzODE(
        const network_type& constants, const state_type& coupling, 
        const state_type& phase_shift, const state_type& omega, unsigned int N, 
        const node_size_type& node_indices)
      :_constants(constants), _coupling(coupling), _phase_shift(phase_shift),
      _omega(omega), _N(N), _node_indices(node_indices) {}


    void operator()(const state_type& x, state_type& dx, const double t) {
      network_type mean_field = CalculateMeanField(x);
      for (size_t i = 0; i < _node_indices.size() - 1; ++i) {
        unsigned int indx = _node_indices[i];
        dx[indx] = (1.-(x[indx]*x[indx]))/2.*mean_field[i][0];
        dx[indx+1] = (1. - x[indx]*x[indx])/(2.*x[indx])*mean_field[i][1];
        dx[indx+2] = (1. + x[indx]*x[indx])/(2.*x[indx])*mean_field[i][1] 
                     + _omega[i];
      }
    }


    network_type CalculateMeanField(const state_type& x) {
      network_type mean_field;
      for (size_t i = 0; i < _node_indices.size() - 1; ++i) {
        mean_field.push_back(state_type(2));
        unsigned int group_indx = _node_indices[i];
        for (size_t j = 0; j < _node_indices.size() - 1; ++j) {
          unsigned int indx = _node_indices[j];
          double cos_factor = cos(_phase_shift[j] - x[group_indx+2]);
          double sin_factor = sin(_phase_shift[j] - x[group_indx+2]);
          double sin_Phi = sin(x[indx+2]), cos_Phi = cos(x[indx+2]);
          double Psi = x[indx+1], rho = x[indx];
          double cos_sum = 0., sin_sum = 0.;
          for (size_t k = 0; k < _constants[j].size(); ++k) {
            double cos_term = (1.+rho*rho)*cos(_constants[j][k]-Psi) + 2.*rho;
            double sin_term = (1.-rho*rho)*sin(_constants[j][k]-Psi);
            double denominator = 1.+rho*rho+2*rho*cos(_constants[j][k]-Psi);
            cos_sum += (cos_term*cos_Phi - sin_term*sin_Phi)/denominator;
            sin_sum += (sin_term*cos_Phi + cos_term*sin_Phi)/denominator;
          }
          mean_field[i][0] += _coupling[j]*(cos_factor*cos_sum 
                                            - sin_factor*sin_sum);
          mean_field[i][1] += _coupling[j]*(sin_factor*cos_sum 
                                            + cos_factor*sin_sum);
        }
        mean_field[i][0] /= static_cast<double>(_N);
        mean_field[i][1] /= static_cast<double>(_N);
      }
      return mean_field;
    }
};
} // nonlinear_systems
#endif
