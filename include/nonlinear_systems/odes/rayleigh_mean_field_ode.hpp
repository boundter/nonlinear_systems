#ifndef __RAYLEIGH_MEAN_FIELD__
#define __RAYLEIGH_MEAN_FIELD__

#include <vector>
#include <stdexcept>

typedef std::vector<double> state_type;

namespace nonlinear_systems {
class RayleighMeanFieldODE {
  public:
    RayleighMeanFieldODE(unsigned int N, const state_type& frequency, 
        double nonlinearity, double coupling)
      :_frequency(frequency), _nonlinearity(nonlinearity), _coupling(coupling),
      _N(N) {
        if (_frequency.size() != _N) {
          throw std::length_error("Length of frequency vector and system size N do not match.");
        } 
    }


  protected:
    state_type _frequency;
    double _nonlinearity, _coupling;
    unsigned int _N;
};

class RayleighMeanFieldODEX: public RayleighMeanFieldODE {
  public:
    RayleighMeanFieldODEX(unsigned int N, const state_type& frequency,
        double nonlinearity, double coupling)
      :RayleighMeanFieldODE(N, frequency, nonlinearity, coupling) {}


    void operator()(const state_type& x, state_type& dx) {
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


class RayleighMeanFieldODEY: public RayleighMeanFieldODE {
  public:
    RayleighMeanFieldODEY(unsigned int N, const state_type& frequency,
        double nonlinearity, double coupling)
      :RayleighMeanFieldODE(N, frequency, nonlinearity, coupling) {}


    void operator()(const state_type& x, state_type& dx) {
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
