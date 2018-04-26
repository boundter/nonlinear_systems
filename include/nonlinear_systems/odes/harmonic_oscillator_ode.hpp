#ifndef __HARMONIC_OSCILLATOR_ODE__
#define __HARMONIC_OSCILLATOR_ODE__

#include <vector>

typedef std::vector<double> state_type;

namespace nonlinear_systems {
/*!
 *  The Harmonic oscillator is described by the ODE
 *  \f[ \ddot{x} = - \omega^2 x. \f]
 *  It has the general analytical solution 
 *  \f[ x(t) = A*\sin(\omega*t + \varphi). \f]
 */
class HarmonicOscillatorODE {
  public:
    double _omega;

    /*!
     *  @param params params will be casted to a double pointer, where the first
     *  entry will be used as the frequency.
     */
    HarmonicOscillatorODE(void* params) {
      _omega = reinterpret_cast<double*>(params)[0];
    }


    void operator()(const state_type& x, state_type& dx, double t) {
      dx[0] = x[1];
      dx[1] = -_omega*_omega*x[0];
    }
};
} // nonlinear_systems

#endif
