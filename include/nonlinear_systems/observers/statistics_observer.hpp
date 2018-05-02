#ifndef __STATISTICS_OBSERVER__
#define __STATISTICS_OBSERVER__

#include <vector>
#include <nonlinear_systems/misc/statistics.hpp>

namespace nonlinear_systems {
template<typename state_type = std::vector<double> >
class AverageObserver {
  public:
    std::vector<double>& _average;
    unsigned int _step_number;

    AverageObserver(std::vector<double>& average)
    : _average(average) {
      _average = std::vector<double>(average.size(), 0);
      _step_number = 1;
    }


    void operator()(const state_type& x, double t) {
      UpdateAverage(x, _step_number, _average);
    }
};


template<typename state_type = std::vector<double> >
class VarianceObserver {
  public:
    std::vector<double>& _average;
    std::vector<double>& _variance;
    unsigned int _step_number;

    VarianceObserver(std::vector<double>& average, std::vector<double>& variance)
    :_average(average), _variance(variance) {
      _step_number = 1;
    }

    
    void operator()(const state_type& x, double t) {
      UpdateSumOfSquares(x, _step_number, _average, _sum_of_squares);
      if (_step_number > 2) {
        for (size_t i = 0; i < _sum_of_squares.size(); ++i) {
          _variance[i] = _sum_of_squares[i]/static_cast<double>(_step_number - 2);
        }
      }
    }

  private:
    std::vector<double> _sum_of_squares;
};
} // nonlinear_systems

#endif
