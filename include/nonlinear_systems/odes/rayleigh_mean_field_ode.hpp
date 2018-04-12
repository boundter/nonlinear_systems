#ifndef __RAYLEIGH_MEAN_FIELD_ODE__
#define __RAYLEIGH_MEAN_FIELD_ODE__

#include <vector>
#include <stdexcept> // std::length_error

typedef std::vector<double> state_type;

namespace nonlinear_systems {
/*!
 * This ODE describes a system Rayleigh oscillators coupled to a mean field. The
 * coupling is the same for all oscillators. *This is not an ODE, but just a
 * template for difference ways to couple to the mean field!*
 */
struct RayleighMeanFieldODE {
    state_type _frequency;
    double _nonlinearity, _coupling;
    unsigned int _N;

    /*!
     *  @param N number of oscillators
     *  @param frequency frequencies for all the oscillators
     *  @param nonlinearity nonlinearity parameter for the system
     *  @param coupling coupling constant of the system
     */
    RayleighMeanFieldODE(unsigned int N, const state_type& frequency, 
        double nonlinearity, double coupling)
      :_frequency(frequency), _nonlinearity(nonlinearity), _coupling(coupling),
      _N(N) {
        if (_frequency.size() != _N) {
          throw std::length_error("Length of frequency vector and system size N do not match.");
        } 
    }


    virtual void operator()(const state_type& x, state_type& dx, const double t)
    {}
};


/*!
 *  This ODE described a system of Rayleigh oscillators coupled to a mean field
 *  in the y-direction. The ODE is 
 *
 *  x'' - nonlinearity*(1-x'^2)*x' + frequency^2*x = coupling*(X - x)
 */
class RayleighMeanFieldODEX: public RayleighMeanFieldODE {
  public:
    
    /*!
     *  @param N number of oscillators
     *  @param frequency frequencies for all the oscillators
     *  @param nonlinearity nonlinearity parameter for the system
     *  @param coupling coupling constant of the system
     */
    RayleighMeanFieldODEX(unsigned int N, const state_type& frequency,
        double nonlinearity, double coupling)
      :RayleighMeanFieldODE(N, frequency, nonlinearity, coupling) {}


    void operator()(const state_type& x, state_type& dx, const double t) {
      double mean_field = CalculateMeanFieldCoordinateX(x);
      for (unsigned int i = 0; i < _N; ++i) {
        dx[2*i] = x[2*i+1];
        dx[2*i+1] = _nonlinearity*(1-x[2*i+1]*x[2*i+1])*x[2*i+1]
          - _frequency[i]*_frequency[i]*x[2*i]
          + _coupling*(mean_field - x[2*i]);
      }
    }


  protected:
  
    double CalculateMeanFieldCoordinateX(const state_type& x) {
      double sum = 0.;
      for (unsigned int i = 0; i < _N; ++i) {
        sum += x[2*i];
      }
      return sum/static_cast<double>(_N);
    }
};


/*!
 *  This ODE described a system of Rayleigh oscillators coupled to a mean field
 *  in the x-direction. The ODE is 
 *
 *  x'' - nonlinearity*(1-x'^2)*x' + frequency^2*x = coupling*(X' - x')
 */
class RayleighMeanFieldODEY: public RayleighMeanFieldODE {
  public:
    
    /*!
     *  @param N number of oscillators
     *  @param frequency frequencies for all the oscillators
     *  @param nonlinearity nonlinearity parameter for the system
     *  @param coupling coupling constant of the system
     */
    RayleighMeanFieldODEY(unsigned int N, const state_type& frequency,
        double nonlinearity, double coupling)
      :RayleighMeanFieldODE(N, frequency, nonlinearity, coupling) {}


    void operator()(const state_type& x, state_type& dx, const double t) {
      double mean_field = CalculateMeanFieldCoordinateY(x);
      for (unsigned int i = 0; i < _N; ++i) {
        dx[2*i] = x[2*i+1];
        dx[2*i+1] = _nonlinearity*(1-x[2*i+1]*x[2*i+1])*x[2*i+1]
          - _frequency[i]*_frequency[i]*x[2*i]
          + _coupling*(mean_field - x[2*i+1]);
      }
    }


  protected:

    double CalculateMeanFieldCoordinateY(const state_type& x) {
      double sum = 0.;
      for (unsigned int i = 0; i < _N; ++i) {
        sum += x[2*i+1];
      }
      return sum/static_cast<double>(_N);
    }
};
}//nonlinear_systems

#endif
