#ifndef __KURAMOTO_SAKAGUCHI_SYSTEM__
#define __KURAMOTO_SAKAGUCHI_SYSTEM__

#include <cmath>
#include <functional>
#include <random>
#include <vector>
#include <nonlinear_systems/systems/generic_system.hpp>
#include <nonlinear_systems/odes/kuramoto_sakaguchi_ode.hpp>

typedef std::vector<double> state_type;

namespace nonlinear_systems {
/*!
 *  \brief Wrapper for the Kuramoto-Sakaguchi ODE.
 *
 *
 *  The system is wrapper for the Kuramoto-Sakaguchi differential equation
 *  \f[ \frac{d\varphi_i}{dt} = \omega_i + \frac{K}{N} \sum_{j=1}^{N} 
 *  \sin(\varphi_j - \varphi_i + \alpha).
 *  \f]
 */
class KuramotoSakaguchiSystem
  : public GenericSystem<KuramotoSakaguchiODE, state_type> {
  public:
    
    /*!
     * @param N number of oscillators.
     * @param frequency the frequency of every single oscillator.
     * @param coupling the coupling between the oscillators.
     * @param phase_shift the phase shift between the oscillators.
     * @param seed the seed for the random number generator.
     */
    KuramotoSakaguchiSystem(unsigned int N, state_type frequency, 
        double coupling, double phase_shift, unsigned long int seed=123456789)
      : GenericSystem<KuramotoSakaguchiODE, state_type>(N, 1) {
      this->_N = N;
      this->_rng.seed(seed);
      this->_ode = std::unique_ptr<KuramotoSakaguchiODE>(
          new KuramotoSakaguchiODE(N, frequency, coupling, phase_shift));
      SetRandomState();
    }

    
    /*!
     *  /brief Distributes the oscillators randomly uniformly along the unit
     *  circle.
     */
    void SetRandomState() {
      std::uniform_real_distribution<double> uniform(-M_PI, M_PI);
      std::function<double()> uniform_dist = std::bind(uniform, 
          std::ref(_rng));
      this->_x = SampleDistribution<state_type, double>(this->_x.size(), 
          &uniform_dist);
    }

    
    /*!
     *  \brief Distributes the oscillators uniformly along the unit circle with
     *  a slight perturbation (one space is left empty).
     */
    void SetPerturbedSplayState() {
      double dist_oscillator = 2*M_PI/(static_cast<double>(N) - 2.);
      for (size_t i = 0; i < this->_x.size(); ++i) {
        this->_x[i] = static_cast<double>(i)*dist_oscillator;
      }
    }

    
    /*!
     *  \brief Distributes the oscillaltors in a cluster with a given width.
     *
     *  @param cluster_width the width of the cluster.
     */
    void SetPerturbedCluster(double cluster_width) {
      double dist_oscillator = cluster_width/(static_cast<double>(N) - 1.);
      for (size_t i = 0; i < this->_x.size(); ++i) {
        this->_x[i] = static_cast<double>(i)*dist_oscillator;
      }
    }

    
    /*!
     *  \brief Calculate the generalized mean field.
     *
     *  The generalized mean field of the order n is defined as 
     *  \f[
     *    Z_n = \frac{1}{N} \sum_{j=1}^{N} e^{i n \varphi_j}.
     *  \f]
     *
     *  @param order the order of the generalized mean field.
     */
    state_type CalculateGeneralizedMeanFieldSpherical(int order) {
      state_type multiple_phases = this->_x;
      for (size_t i = 0; i < multiple_phases.size(); ++i) {
        multiple_phases[i] = static_cast<double>(order)*(multiple_phases % (2*M_PI));
      }
      return CalculateMeanFieldSpherical(multiple_phases.begin(), 
          multiple_phases.end());
    }

  protected:
    std::mt19937_64 _rng;
}
} // nonlinear_systems

#endif
