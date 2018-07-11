#ifndef __REDUCED_M_KURAMOTO_SAKAGUCHI_WATANABE_STROGATZ_ODE__
#define __REDUCED_M_KURAMOTO_SAKAGUCHI_WATANABE_STROGATZ_ODE__

#include <cmath> // sin,cos
#include <vector>

typedef std::vector<double> state_type;
typedef std::vector<state_type> network_type;

namespace nonlinear_systems {

class ReducedMKuramotoSakaguchiWatanabeStrogatzODE {
  public:
    network_type _constants;
    state_type _coupling, _phase_shift;
    unsigned int _N;


    ReducedMKuramotoSakaguchiWatanabeStrogtazODE(
        const network_type& constants, const state_type& coupling, 
        const state_type& phase_shift, unsigned int N)
      :_constants(constants), _coupling(coupling), phase_shift(phase_shift),
      _N(N) {}


    void operator()(const state_type& x, state_type& dx, const double t) {
      network_type mean_field = CalculateMeanField(x);
      // TODO: Check for correct size
      for (size_t i = 0; i < x.size()/3; ++i) {
        dx[i*3] = (1.-(x[i*3]*x[i*3]))/2.*mean_field[i][0];
        dx[i*3+1] = (1. - x[i*3]*x[i*3])/(2.*x[i*3])*mean_field[i][1];
        dx[i*3+2] = (1. + x[i*3]*x[i*3])/(2.*x[i*3])*mean_field[i][1] + omega[i];
      }
    }


    network_type CalculateMeanField(const state_type& x) {
      network_type mean_field;
      for (size_t i = 0; i < x.size()/3, ++i) {
        mean_field.push_back(state_type(2));
        for (size_t j = 0; j < x.size()/3; ++j) {
          double cos_factor = cos(_phase_shift[j] - x[i*3+2]);
          double sin_factor = sin(_phase_shift[j] - x[i*3+2]);
          double sin_Phi = sin(x[3*j+2]), cos_Phi = cos(x[3*j+2]);
          double Psi = x[3*j+1], rho = x[3*j];
          double cos_sum = 0., sin_sum = 0.;
          for (size_t k = 0; k < _constants[j].size(); ++k) {
            double cos_term = (1.+rho*rho)*cos(_constants[i][j]-Psi) + 2.*rho;
            double sin_term = (1.-rho*rho)*sin(_constants[i][j]-Psi);
            double denominator = 1.+rho*rho+2*rho*cos(_constants[i][j]-Psi);
            cos_sum += (cos_term*cos_Phi - sin_term*sin_Phi)/denominator;
            sin_sum += (sin_term*cos_Phi + cos_term*sin_Phi)/denominator;
          }
          mean_field[i][0] += _coupling[j]*(cos_factor*cos_sum 
                                            - sin_factor*sin_sum);
          mean_field[i][1] += _coupling[j]*(sin_factor*cos_sum 
                                            + cos_factor*sin_sum);
        }
        mean_field[i][0] /= _N;
        mean_field[i][1] /= _N;
      }
      return mean_field;
    }
};
} // nonlinear_systems
#endif
