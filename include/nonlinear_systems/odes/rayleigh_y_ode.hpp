#ifndef __RAYLEIGH_Y_ODE__
#define __RAYLEIGH_Y_ODE__

#include <vector>

namespace nonlinear_systems {
  namespace ode {
    template<typename state_type = std::vector<double> >
      class RayleighYODE {
        public:
          state_type omega;
          unsigned int N;
          double eta, eps;


          RayleighYODE(const state_type& frequency, unsigned int N_oscillators, 
              double nonlinearity, double coupling)
            : omega(frequency), N(n_oscillator), eta(nonlinearity), 
            eps(coupling) { };


          double CalculateMeanFieldY(const state_type& y) {
            double sum = 0.;
            for (unsigned int i = 0; i < N; ++i) {
              sum += y[2*i+1];
            }
            return sum/(double)N;
          };


          void operator() (const state_type& y, state_type& dy, 
              const double t) {
            double Y = CalculateMeanFieldY(y);
            for (unsigned int i = 0; i < N; ++i) {
              dy[2*i] = y[2*i+1];
              dy[2*i+1] = eta*(1-y[2*i+1]*y[2*i+1])*y[2*i+1] 
                - omega[i]*omega[i]*y[2*i] + eps*(Y-y[2*i+1]);
            }
          };
      };
  } // ode
} // nonlinear_systems

#endif
