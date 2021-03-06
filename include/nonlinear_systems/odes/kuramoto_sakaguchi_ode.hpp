#ifndef __KURAMOTO_SAKAGUCHI_ODE__
#define __KURAMOTO_SAKAGUCHI_ODE__

#include <cmath>
#include <stdexcept>
#include <vector>

typedef std::vector<double> state_type;

namespace nonlinear_systems {
/*!
 *  ODE for the Kuramoto-Sakaguchi model. It reads
 * /f[ \dot{\varphi}_i = \omega_i + \frac{K}{N} \sum_{j=1}^{N} \sin(\varphi_j -
 * \varphi_i + \alpha). /f]
 */
class KuramotoSakaguchiODE {
  public:
    unsigned int _N;
    state_type _frequency;
    double _coupling;
    double _phase_shift;

    /*!
     *  @param N number of oscillators.
     *  @param frequency the natural frequency of the oscillators.
     *  @param coupling the coupling constant.
     *  @param phase_shift the phase shift of the ode.
     */
    KuramotoSakaguchiODE(unsigned int N, state_type frequency, double coupling, 
        double phase_shift)
   :_N(N), _frequency(frequency), _coupling(coupling), _phase_shift(phase_shift) {
      if (_frequency.size() != _N) {
        throw std::length_error("Length of frequency vector and system size N do not match.");
      }
    }

    // solves the integration with the mean field the new ODE reads
    // \dot{\varphi}_i = \omega_i + R*K\sin(\Psi - \varphi_i + \alpha)
    void operator()(const state_type& x, state_type& dx, double t) {
      state_type mean_field = CalculateMeanField(x);
      for (unsigned int i = 0; i < _N; ++i) {
        dx[i] = _frequency[i] + _coupling*mean_field[0]
                                *sin(mean_field[1] - x[i] + _phase_shift);
      }
    }

    /*!
     *  Returns the mean field \f$ Z = R e^{i\Psi} \f$ as a vector of the form
     *  {R, Psi}. 
     */
    state_type CalculateMeanField(const state_type& x) {
      state_type mean_field(2);
      double re = 0., im = 0.;
      for (unsigned int i = 0; i < _N; ++i) {
        re += cos(x[i]);
        im += sin(x[i]);
      }
      re /= static_cast<double>(_N);
      im /= static_cast<double>(_N);
      mean_field[0] = sqrt(re*re + im*im);
      mean_field[1] = atan2(im, re);
      return mean_field;
    }
};
} // nonlinear_systems

#endif
