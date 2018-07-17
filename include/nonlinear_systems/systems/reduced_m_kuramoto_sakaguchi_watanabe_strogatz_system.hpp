#ifndef __REDUCED_M_KURAMOTO_SAKAGUCHI_WATANABE_STROGATZ_SYSTEM__
#define __REDUCED_M_KURAMOTO_SAKAGUCHI_WATANABE_STROGATZ_SYSTEM__

#include <cmath>
#include <exception>
#include <iostream>
#include <vector>
#include <nonlinear_systems/odes/reduced_m_kuramoto_sakaguchi_watanabe_strogatz_ode.hpp>
#include <nonlinear_systems/systems/generic_network.hpp>

#include <nlopt.hpp>

typedef std::vector<double> state_type;
typedef std::vector<state_type> network_type;
typedef std::vector<unsigned int> node_size_type;

namespace nonlinear_systems {
class ReducedMKuramotoSakaguchiWatanabeStrogatzSystem
  : public GenericNetwork<ReducedMKuramotoSakaguchiWatanabeStrogatzODE, double> {
  public:

    ReducedMKuramotoSakaguchiWatanabeStrogatzSystem(
        const network_type& phases, const char* conversion = "splay")
   :GenericNetwork<ReducedMKuramotoSakaguchiWatanabeStrogatzODE, double>(
       node_size_type(phases.size(), 1), 3) {
      this->_x.resize(3*phases.size());
      if (conversion == "splay") {
        TransformPhasesToWS(phases);
      }
      else if (conversion == "identity") {
        IdentityConversion(phases);
      }
      else {
        throw std::invalid_argument("Unknown conversion.");
      }
      /*
      this->_ode = std::unique_ptr<ReducedMKuramotoSakaguchiWatanabeStrogatzODE>(
          new ReducedMKuramotoSakaguchiWatanabeStrogatzODE());
      */
    }


    network_type CalculatePhases() {
      network_type phases;
      for (size_t i = 0; i < _constants.size(); ++i) {
        double rho = this->_x[_node_indices[i]];
        double Psi = this->_x[_node_indices[i]+1];
        double Phi = this->_x[_node_indices[i]+2];
        phases.push_back(TransformConstantsToPhases(_constants[i], rho, Psi, 
              Phi));
      }
      return phases;
    }
  
  
  protected:
    network_type _constants;


    void IdentityConversion(const network_type& phases) {
      for (size_t i = 0; i < phases.size(); ++i) {
        this->_x[_node_indices[i]] = 0.;
        this->_x[_node_indices[i]+1] = 0.;
        this->_x[_node_indices[i]+2] = 0.;
        _constants.push_back(phases[i]);
      }
    }


    state_type TransformConstantsToPhases(const state_type& constants,
        double rho, double Psi, double Phi) {
      state_type phases(constants.size());
      for (size_t i = 0; i < constants.size(); ++i) {
        double cos_term = (1+rho*rho)*cos(constants[i]-Psi) + 2*rho;
        double sin_term = (1-rho*rho)*sin(constants[i]-Psi);
        double im = sin_term*cos(Phi) + cos_term*sin(Phi);
        double re = cos_term*cos(Phi) - sin_term*sin(Phi);
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
   
    
    // x[0] = rho, x[1] = Phi
    static double CalculatePotential(const state_type& x, state_type& grad, 
        void* args) {
      // Does this work?
      state_type phases = *reinterpret_cast<state_type*>(args);
      double sum = 0.;
      for (size_t i = 0; i < phases.size(); ++i) {
        sum += log((1 + x[0]*x[0] - 2*x[0]*cos(phases[i]-x[1]))
                   /((1. - x[0]*x[0])/2.));
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
      for (size_t i = 0; i < phases.size(); ++i) {
        // Initial guess
        state_type mean_field = CalculateKuramotoMeanField(phases[i]);
        // minimize the potential
        nlopt::opt opt(nlopt::LN_SBPLX, 2);
        state_type lower_bounds = {0., -M_PI};
        state_type upper_bounds = {0.99, M_PI};
        opt.set_lower_bounds(lower_bounds);
        opt.set_upper_bounds(upper_bounds);
        state_type current_phases = phases[i];
        opt.set_min_objective(CalculatePotential, &current_phases);
        opt.set_xtol_rel(1e-6);
        double minimum;
        try {
          nlopt::result result = opt.optimize(mean_field, minimum);
        }
        catch(std::exception &e) {
          std::cout << "Could not find minimum of potential for group " << i <<
            " with nlopt error: " << e.what() << std::endl;
        }
        double rho = mean_field[0];
        double Phi = mean_field[1];
        double Psi = CalculatePsi(current_phases, rho, Phi);
        // Make use of the transformation phases<->constants, Phi<->Psi,
        // rho->-rho, which changes the role of phases and constants
        _constants.push_back(TransformConstantsToPhases(current_phases, -rho, Phi, 
              Psi));
        this->_x[_node_indices[i]] = rho;
        this->_x[_node_indices[i]+1] = Psi;
        this->_x[_node_indices[i]+2] = Phi;
      }
    }
};
} // nonlinear_systems
#endif
