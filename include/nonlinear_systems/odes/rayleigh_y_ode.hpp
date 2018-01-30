#ifndef __RAYLEIGH_Y_ODE__
#define __RAYLEIGH_Y_ODE__

#include <vector>

namespace nonlinear_systems {
  template<typename state_type = std::vector<double> >
    class RayleighYODE {
      public:
        state_type omega;
        unsigned int N;
        double eta, eps;

        RayleighYODE(const state_type& frequency, unsigned int N_oscillators, 
            double nonlinearity, double coupling);
        double CalculateMeanFieldY(const state_type& y);
        void operator() (const state_type& y, state_type& dy, 
            const double t);
    };
} // nonlinear_systems

using namespace nonlinear_systems;

template<typename state_type>
RayleighYODE<state_type>::RayleighYODE(const state_type& frequency, 
    unsigned int N_oscillator, double nonlinearity, double coupling)
  : omega(frequency), N(N_oscillator), eta(nonlinearity), eps(coupling) {}

template<typename state_type>
double RayleighYODE<state_type>::CalculateMeanField(const state_type& x) {
  double sum = 0.;
  for (unsigned int i = 0; i < N; ++i) {
    sum += x[2*i+1];
  }
  return sum/static_cast<double>(N);
}

template<typename state_type>
void RayleighYODE<state_type>::operator()(const state_type& x, state_type& dx,
    const double t) {
  double Y = CalculateMeanFieldY(y);
  for (unsigned int i = 0; i < N; ++i) {
    dx[2*i] = x[2*i+1];
    dx[2*i+1] = eta*(1-x[2*i+1]*x[2*i+1])*x[2*i+1] 
              - omega[i]*omega[i]*x[2*i] + eps*(Y-x[2*i+1]);
  } 
}
#endif
