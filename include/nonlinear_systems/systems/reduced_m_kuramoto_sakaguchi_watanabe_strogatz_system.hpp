#ifndef __REDUCED_M_KURAMOTO_SAKAGUCHI_WATANABE_STROGATZ_SYSTEM__
#define __REDUCED_M_KURAMOTO_SAKAGUCHI_WATANABE_STROGATZ_SYSTEM__

#include <cmath>
#include <vector>
#include <nonlinear_systems/odes/reduced_m_kuramoto_sakaguchi_watanabe_strogatz_ode.hpp>
#include <nonlinear_systems/systems/generic_network.hpp>

typedef std::vector<double> state_type;
typedef std::vector<state_type> network_type;
typedef std::vector<unsigned int> node_size_type;

namespace nonlinear_systems {
class ReducedMKuramotoSakaguchiWatanabeStrogatzSystem
  : public GenericNetwork<ReducedMKuramotoSakaguchiWatanabeStrogatzODE, double> {
  public:

    ReducedMKuramotoSakaguchiWatanabeStrogatzSystem(
        const network_type& phases)
   :GenericNetwork<ReducedMKuramotoSakaguchiWatanabeStrogatzODE, double>(
       node_size_type(phases.size()), 3) {
      TransformPhasesToWS(phases);
      /*
      this->_ode = std::unique_ptr<ReducedMKuramotoSakaguchiWatanabeStrogatzODE>(
          new ReducedMKuramotoSakaguchiWatanabeStrogatzODE());
      */
    }


    network_type CalculatePhases() {
      network_type phases(_constants.size());
      for (size_t i = 0; i < _constants.size(); ++i) {
        phases[i].push_back(state_type(_constants[i].size()));
        double rho = x[_node_indices[i]];
        double Psi = x[_node_indices[i+1]];
        double Phi = x[_node_indices[i+2]];
        phases[i].push_back(TransformConstantsToPhases(_constants[i], rho, Psi, 
              Phi));
      }
      return phases;
    }
  
  
  protected:
    network_type _constants;


    state_type TransformConstantsToPhases(const state_type& constants,
        double rho, double Psi, double Phi) {
      state_type phases(constants.size());
      for (size_t i = 0; i < constants.size(); ++i) {
        double cos_term = (1+rho*rho)*cos(constants[i]-Psi) + 2*rho;
        double sin_term = (1-rho*rho)*sin(constants[i]-Psi);
        double im = sin_term*cos(Phi) + cos_term*sin(Phi);
        double re = cos_term*cos(Phi) - sin_term*cos(Phi);
        phases[i] = atan2(im, re);
      }
      return phases;
    }

    state_type CalculateKuramotoMeanField(const state_type& phases) {
      double re = 0., im = 0.;
      for (size_t i = 0; i < phases.size(); ++i) {
        re += cos(phases[i]);
        im += sin(phases[i]);
      }
      state_type mean_field(2);
      mean_field[0] = sqrt(re*re + im*im)/static_cast<double>(phases.size());
      mean_field[1] = atan2(im, re);
      return mean_field;
    }
   

    double CalculatePotential(const state_type& phases, double rho, double Phi) {
      double sum = 0.;
      for (size_t i = 0; i < phases.size(); ++i) {
        sum += log((1. + rho*rho - 2.*rho*cos(phases[i]-Phi))
                   /((1. - rho*rho)/2.));
      }
      return sum/static_cast<double>(phases.size());
    }


    double CalculatePsi(const state_type& phases, double rho, double Phi) {
      double sum = 0.;
      for (size_t i = 0; i < phases.size(); ++i) {
        sum += atan((1.+rho)/(1.-rho)*tan((phases[i]-Phi)/2.));
      }
      return -2.*sum/static_cast<double>(phases.size());
    }


    void TransformPhasesToWS(const network_type& phases) {
      _constants.resize(phases.size());
      for (size_t i = 0; i < phases.size(); ++i) {
        state_type mean_field = CalculateKuramotoMeanField(phases[i]);
        // Initial guess
        double rho = mean_field[0];
        double Phi = mean_field[1];
        // TODO: Somehow the potential is minimizmed and we get the correct rho, Phi
        double Psi = CalculatePsi(phases, rho, Phi);
        // Make use of the transformation phases<->constants, Phi<->Psi,
        // rho->-rho, which changes the role of phases and constants
        constants[i].push_back(TransformConstantsToPhases(phases, -rho, Phi, 
              Psi));
        this->_x[this->_node_indices[i]] = rho;
        this->_x[this->_node_indices[i]+1] = Psi;
        this->_x[this->_node_indices[i]+2] = Phi;
      }
    }
};
} // nonlinear_systems
#endif
