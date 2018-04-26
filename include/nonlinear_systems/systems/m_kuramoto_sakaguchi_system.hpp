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
/*!
 * The M-Kuramoto-Sakaguchi ODE describes M groups of Kuramoto-Sakaguchi
 * oscillators, where each group can have a different coupling and phase shift.
 * The ODE for an oscillator reads
 * \f[ \dot{\varphi}_i^{\sigma} = \omega_i + \sum_{\sigma' = 1}^M
 * \frac{K_{\sigma\sigma'}}{N} \sum_{j=1}^{N_{\sigma'}} 
 * \sin(\varphi_j^{\sigma'} - \varphi_i^{\sigma} + \alpha_{\sigma\sigma'}) \f]
 */
class MKuramotoSakaguchiSystem
  : public GenericNetwork<MKuramotoSakaguchiODE, double> {
  public:

    /*!
     * @param frequency flattened representation of the oscillators natural
     * frequency.
     * @param coupling matrix of the coupling between groups.
     * @param phase_shift matrix of the phase shift between groups.
     * @param node_size the size of the groups/nodes.
     * @param seed the seed used for the internal random number generator.
     */
    MKuramotoSakaguchiSystem(const state_type& frequency, 
        const network_type& coupling, const network_type& phase_shift,
        const node_size_type& node_size, unsigned int seed=123456789)
      :GenericNetwork<MKuramotoSakaguchiODE, double>(node_size, 1) {
      _rng.seed(seed);
      SetRandomUniformState();
      _coupling = coupling;
      _phase_shift = phase_shift;
      _node_size = node_size;
      this->_ode = std::unique_ptr<MKuramotoSakaguchiODE>(
          new MKuramotoSakaguchiODE(frequency, _coupling, _phase_shift, 
            this->_node_indices));
    }


    /*!
     * A simplified constructor for a system of two groups of identical
     * oscillators without phase shift, where the only depends on the acting
     * group, so \f$ K_{\sigma\sigma'} = K_{\sigma'} \f$. The system can be
     * rotated, so that the first group has no natural frequency, and the time
     * scaled so that the first group has a coupling strength of 1. This way the
     * system can be parametrized by two numbers.
     *
     * @param frequency the natural frequency of the second group.
     * @param repulsive_excess the repulsive excess \f$ \varepsilon = -(1+K_2) \f$.
     * @param node_size the size of the groups/nodes.
     * @param seed the seed for the internal random number generator.
     */
    MKuramotoSakaguchiSystem(double frequency, double repulsive_excess,
        const node_size_type& node_size, unsigned int seed=123456789)
      :GenericNetwork<MKuramotoSakaguchiODE, double>(node_size, 1) {
      _rng.seed(seed); 
      SetRandomUniformState();
      state_type frequency_vector(this->_node_indices.back());
      for (size_t i = this->_node_indices[1]; i < frequency_vector.size(); ++i) {
        frequency_vector[i] = frequency;   
      }
      // defining this also implicitly checks for the correct system size
      state_type coupling_row = {1., -(1.+repulsive_excess)};
      network_type coupling = {coupling_row, coupling_row};
      state_type zero_row = {0., 0.};
      network_type phase_shift = {zero_row, zero_row};
      this->_ode = std::unique_ptr<MKuramotoSakaguchiODE>(
          new MKuramotoSakaguchiODE(frequency_vector, coupling, phase_shift, 
            this->_node_indices));
      
    }


    /*!
     *  Distributes the phases uniform randomly along the unit circle.
     */
    void SetRandomUniformState() {
      std::uniform_real_distribution<double> uniform(-M_PI, M_PI);
      std::function<double()> uniform_dist = std::bind(uniform, std::ref(_rng));
      this->_x = SampleDistribution<state_type, double>(this->_x.size(), 
          &uniform_dist);
    }


    /*!
     *  Distributes the phases around clusters. Careful: If
     *  clustre_distance*number_clusters > 2*pi the clusters will be wrapped
     *  around!
     *
     *  @param cluster_width the width of every single cluster, oscillators will
     *  be placed with uniform distance to eachother.
     *  @param cluster_distance the distance between clusters.
     */
    void SetPerturbedClusters(double cluster_width, double cluster_distance) {
      size_t number_clusters = _node_size.size();
      state_type cluster_position;
      for (size_t i = 0; i < number_clusters; ++i) {
        cluster_position.push_back(i*cluster_distance);
      }
      for (size_t i = 0; i < this->_node_indices.size() - 1; ++i) {
        double distance_between_oscillators = 
          cluster_width/(static_cast<double>(_node_size[i]) - 1.);
        for (unsigned int j = this->_node_indices[i]; 
            j < this->_node_indices[i+1]; ++j) {
          this->_x[j] = cluster_position[i] - cluster_width/2. 
            + distance_between_oscillators*(j - this->_node_indices[i]);
        }
      }
    }


    /*!
     *  Calculates the mean field for the seperate groups. The first index gives
     *  the group number and the second one the quantity; 0 is the order
     *  parameter and 1 the phase.
     */
    network_type CalculateMeanField() {
      return this->_ode->CalculateMeanField(this->_x);
    }


    /*!
     *  Calculates the generalized mean fields/ fourier coefficients. The format
     *  is like for the mean field.
     *
     *  @param fourier_number the number of the generlaized mean field to
     *  calculate.
     */
    network_type CalculateGeneralizedMeanField(int fourier_number) {
      state_type new_state(this->_x.size());
      for (size_t i = 0; i < new_state.size(); ++i) {
        new_state[i] = static_cast<double>(fourier_number)*this->_x[i];
      }
      return this->_ode->CalculateMeanField(new_state);
    }


    // TODO: This should be 3-dimensional, \sigma, \sigma', complex number
    /*!
     *  The forcing generalizes the mean field for the M-Kuramoto-Sakaguchi
     *  system. With it the differential equation can be reduced. The forcing is 
     *  \f[ H_{\sigma\sigma'} = \sum_{\sigma'} K_{\sigma\sigma'} 
     *  \frac{N_\sigma'}{N} Z_{\sigma'}e^{i\alpha_{\sigma\sigma'}}. \f]
     *
     *  This seems to be wrong right now.
     */
    network_type CalculateForcing() {
      network_type forcing(_node_size.size());
      network_type mean_field = CalculateMeanField();
      for (size_t i = 0; i < forcing.size(); ++i) {
        double real_part = 0., imag_part = 0.;
        for (size_t j = 0; j < _node_size.size(); ++j) {
          double coupling_part = _coupling[i][j]*_node_size[j]
            /this->_node_indices.back()*mean_field[j][0];
          real_part += coupling_part*cos(mean_field[j][1]+_phase_shift[i][j]);
          imag_part += coupling_part*sin(mean_field[j][1]+_phase_shift[i][j]);
        }
        forcing[i].push_back(sqrt(real_part*real_part + imag_part*imag_part));
        forcing[i].push_back(atan2(imag_part, real_part));
      }
      return forcing;
    }


  protected:
    std::mt19937_64 _rng;
    network_type _coupling, _phase_shift;
    node_size_type _node_size;
    
};
} // nonlinear_systems
#endif
